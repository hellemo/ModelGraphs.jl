#####################################################
# Link Model
# A link model stores model data about graph links
#####################################################
mutable struct LinkModel <: JuMP.AbstractModel   #subtyping here so I can get ConstraintRef

    graphvariables::Dict{Int,GraphVariable}              #global level variables.  Defined over an entire graph
    graphconstraints::Dict{Int,JuMP.AbstractConstraint}  #global constraint.  Defined over graph variables.  Think "first stage" constraints
    linkconstraints::Dict{Int,LinkConstraint}            #linking constraint.  Defined over variables in nodes

    graphvarnames::Dict{Int,String}
    graphconstraintnames::Dict{Int,String}
    linkconstraintnames::Dict{Int,String}

    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar
    solver::Union{MOI.AbstractOptimizer,AbstractGraphSolver}
end

#Constructor
#LinkModel(;solver = MOI.) = LinkModel(LinkData(),0,JuMP.AffExpr(),solver)

#Need this to dispatch @constraint with graph variables
struct GraphVariable <: JuMP.AbstractVariable
end

#Reference to a Graph Variable.  These variables can be referenced in constraints from JuMP models
struct GraphVariableRef <: JuMP.AbstractVariableRef
    model::LinkModel
    idx::Int
end

#TODO: I think I need a custom variable type to dispatch on @constraint for node models.  Node models could add constraints that use the graph model.
Base.broadcastable(v::GraphVariableRef) = Ref(v)
Base.copy(v::GraphVariableRef) = v
Base.:(==)(v::GraphVariableRef, w::GraphVariableRef) = v.model === w.model && v.idx == w.idx
JuMP.owner_model(v::GraphVariableRef) = v.model
JuMP.isequal_canonical(v::GraphVariableRef, w::GraphVariableRef) = v == w
JuMP.variable_type(::LinkModel) = GraphVariableRef
function JuMP.add_variable(m::LinkModel, v::JuMP.AbstractVariable, name::String="")
    m.nextvaridx += 1
    vref = GraphVariableRef(m, m.nextvaridx)
    m.variables[vref.idx] = v
    JuMP.set_name(vref, name)
    vref
end

function MOI.delete!(m::LinkModel, vref::GraphVariableRef)
    delete!(m.linkvariables, vref.idx)
    delete!(m.linkvarnames, vref.idx)
end
MOI.is_valid(m::LinkModel, vref::GraphVariableRef) = vref.idx in keys(m.linkvariables)
JuMP.num_variables(m::LinkModel) = length(m.linkvariables)

#####################################################
#   Link Constraint
#   A linear constraint between JuMP Models (nodes)
#   Link constraints can be equality or inequality.
#####################################################
abstract type LinkSet <: MOI.AbstractScalarSet end

struct LinkConstraint{F <: AbstractJuMPScalar,S <: LinkSet} <: AbstractConstraint
    func::F
    set::S
end

function JuMP.build_constraint(_error::Function, func::AbstractJuMPScalar,set::LinkSet)
    constraint = LinkConstraint(func, set)
    return constraint
end

jump_function(constraint::LinkConstraint) = constraint.func
moi_set(constraint::LinkConstraint) = constraint.set
shape(::LinkConstraint) = JuMP.ScalarShape()

#TODO Figuring out how this will work
function constraint_object(ref::ConstraintRef{Model, _MOICON{FuncType, LinkSet}}) where
        {FuncType <: MOI.AbstractScalarFunction}
    model = ref.model
    f = MOI.get(model, MOI.ConstraintFunction(), ref)::FuncType
    s = MOI.get(model, MOI.ConstraintSet(), ref)::LinkSet
    return ScalarConstraint(jump_function(model, f), s)
end

# Constraints
struct GraphConstraintRef
    model::LinkModel # `model` owning the constraint
    idx::Int         # Index in `model.constraints`
end
JuMP.constraint_type(::LinkModel) = GraphConstraintRef

#Add a GraphConstraint to a LinkModel
function JuMP.add_constraint(m::LinkModel, con::JuMP.AbstractConstraint, name::String="")
    index = getnumgraphconstraints(m) + 1
    cref = GraphConstraintRef(m, index)
    m.graphconstraints[cref.idx] = con
    JuMP.set_name(cref, name)
    return cref
end

#Add a LinkConstraint to a LinkModel
function JuMP.add_constraint(m::LinkModel, con::LinkConstraint, name::String="")
    index = getnumlinkconstraints(m) + 1
    cref = GraphConstraintRef(m, index)
    m.linkconstraints[cref.idx] = con
    JuMP.set_name(cref, name)
    return cref
end

function MOI.delete!(m::LinkModel, cref::GraphConstraintRef)
    if typeof(constraint_object(cref)) == LinkConstraint
        delete!(m.linkconstraints, cref.idx)
        delete!(m.linkconnames, cref.idx)
    else
        delete!(m.graphconstraints, cref.idx)
        delete!(m.graphconnames, cref.idx)
    end
end

MOI.is_valid(m::LinkModel, cref::GraphConstraintRef) = cref.idx in keys(m.linkconstraints) || cref.idx in keys(m.graphconstraints)
function JuMP.constraint_object(cref::GraphConstraintRef, F::Type, S::Type)
    con = cref.model.constraints[cref.idx]
    con.func::F
    con.set::S
    return con
end

#Get the  nodes in a link constraint
function StructureGraphs.getnodes(con::LinkConstraint)
    vars = keys(con.func.terms)
    nodes = unique([getnode(var) for var in vars])  #unique might lead to slow behavior
    return nodes
end

function getnumnodes(con::LinkConstraint)
    nodes = getnodes(con)
    return length(nodes)
end


#getlinkdata(model::LinkModel) = model.linkdata

#Get the 2 variable or multi-variable linkconstraints
#getlinkconstraints(model::LinkModel) = getlinkdata(model).linkconstraints
#getsimplelinkconstraints(model::LinkModel) = getlinkdata(model).linkconstraints[getlinkdata(model).simple_links]
#gethyperlinkconstraints(model::LinkModel) = getlinkdata(model).linkconstraints[getlinkdata(model).hyper_links]

is_simplelinkconstr(con::LinkConstraint) = getnumnodes(con) == 2 ? true : false
is_hyperlinkconstr(con::LinkConstraint) = getnumnodes(con) > 2 ? true : false

#Extend JuMP's add constraint for link models.  Return a reference to the constraint
#NOTE JuMP no longer uses LinearConstraint
function JuMP.addconstraint(model::LinkModel,constr::JuMP.LinearConstraint)
    #TODO Do some error checking here
    linkdata = getlinkdata(model)
    linkconstr = LinkConstraint(constr)
    push!(linkdata.linkconstraints,linkconstr)
    ref = ConstraintRef{LinkModel,LinkConstraint}(model, length(linkdata.linkconstraints))

    if is_linkconstr(linkconstr)
        push!(linkdata.simple_links,length(linkdata.linkconstraints))
    elseif is_hyperconstr(linkconstr )
        push!(linkdata.hyper_links,length(linkdata.linkconstraints))
    else
        error("constraint $constr doesn't make sense")
    end
    return ref
end


#Constructor
#NOTE Likely not used anymore
# LinkConstraint(con::JuMP.LinearConstraint) = LinkConstraint(con.terms,con.lb,con.ub)
# LinkConstraint(ref::ConstraintRef) = ref.m.linkdata.linkconstraints[ref.idx]  #Get the Link constraint from a constraint reference
# LinkConstraint(con::LinkConstraint) = LinkConstraint(con.terms,con.lb,con.ub)
