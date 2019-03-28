import MathProgBase.SolverInterface:AbstractMathProgSolver

struct LinkVariableRef <: JuMP.AbstractVariableRef
    model::LinkModel
    idx::Int
end

Base.broadcastable(v::LinkVariableRef) = Ref(v)
Base.copy(v::LinkVariableRef) = v
Base.:(==)(v::LinkVariableRef, w::LinkVariableRef) = v.model === w.model && v.idx == w.idx
JuMP.owner_model(v::LinkVariableRef) = v.model
JuMP.isequal_canonical(v::LinkVariableRef, w::LinkVariableRef) = v == w
JuMP.variable_type(::LinkModel) = LinkVariableRef
function JuMP.add_variable(m::LinkModel, v::JuMP.AbstractVariable, name::String="")
    m.nextvaridx += 1
    vref = LinkVariableRef(m, m.nextvaridx)
    m.variables[vref.idx] = v
    JuMP.set_name(vref, name)
    vref
end

function MOI.delete!(m::LinkModel, vref::LinkVariableRef)
    delete!(m.linkvariables, vref.idx)
    delete!(m.linkvarnames, vref.idx)
end
MOI.is_valid(m::LinkModel, vref::LinkVariableRef) = vref.idx in keys(m.linkvariables)
JuMP.num_variables(m::StructuredModel) = length(m.variables)

#####################################################
#   Link Constraint
#   A linear constraint between JuMP Models (nodes)
#   Link constraints can be equality or inequality.
#####################################################
# mutable struct LinkConstraint <: JuMP.AbstractConstraint
#     terms::JuMP.AffExpr
#     lb::Number
#     ub::Number
# end

#Constructor
# LinkConstraint(con::JuMP.LinearConstraint) = LinkConstraint(con.terms,con.lb,con.ub)
# LinkConstraint(ref::ConstraintRef) = ref.m.linkdata.linkconstraints[ref.idx]  #Get the Link constraint from a constraint reference
# LinkConstraint(con::LinkConstraint) = LinkConstraint(con.terms,con.lb,con.ub)

# #####################################################
# # Link Data
# #####################################################
# mutable struct LinkData
#     linkconstraints::Vector{LinkConstraint}          #all links
#     simple_links::Vector{Int}  #references to the 2 node link constraints
#     hyper_links::Vector{Int}   #references to linkconstraints with 3 or more nodes
# end
# LinkData() = LinkData(Vector{LinkConstraint}(),Vector{Int}(),Vector{Int}())

# Constraints
struct LinkConstraintRef
    model::StructuredModel # `model` owning the constraint
    idx::Int       # Index in `model.constraints`
end
JuMP.constraint_type(::LinkModel) = LinkConstraintRef
function JuMP.add_constraint(m::LinkModel, c::JuMP.AbstractConstraint, name::String="")
    m.nextconidx += 1
    cref = LinkConstraintRef(m, m.nextconidx)
    m.linkconstraints[cref.idx] = c
    JuMP.set_name(cref, name)
    cref
end
function MOI.delete!(m::LinkModel, cref::LinkConstraintRef)
    delete!(m.constraints, cref.idx)
    delete!(m.connames, cref.idx)
end
MOI.is_valid(m::LinkModel, cref::LinkConstraintRef) = cref.idx in keys(m.constraints)
function JuMP.constraint_object(cref::LinkConstraintRef, F::Type, S::Type)
    c = cref.model.constraints[cref.idx]
    c.func::F
    c.set::S
    c
end

#Get the  nodes in a link constraint
function StructureGraphs.getnodes(con::LinkConstraint)
    vars = con.terms.vars
    nodes = unique([getnode(var) for var in vars])
    return nodes
end

function getnumnodes(con::LinkConstraint)
    nodes = getnodes(con)
    return length(nodes)
end

#####################################################
# Link Model
#   A link model stores model data about graph links
#####################################################
mutable struct LinkModel <: JuMP.AbstractModel   #subtyping here so I can get ConstraintRef
    linkdata::LinkData           #LinkModel's store indices for each link constraint added

    linkconstraints
    simplelinks
    hyperlinks

    linkvariables
    linkvarnames

    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar
    solver::Union{AbstractMathProgSolver,AbstractGraphSolver}
end
LinkModel(;solver = JuMP.UnsetSolver()) = LinkModel(LinkData(),0,JuMP.AffExpr(),solver)
getlinkdata(model::LinkModel) = model.linkdata

#Get the 2 variable or multi-variable linkconstraints
getlinkconstraints(model::LinkModel) = getlinkdata(model).linkconstraints
getsimplelinkconstraints(model::LinkModel) = getlinkdata(model).linkconstraints[getlinkdata(model).simple_links]
gethyperlinkconstraints(model::LinkModel) = getlinkdata(model).linkconstraints[getlinkdata(model).hyper_links]

is_simplelinkconstr(con::LinkConstraint) = getnumnodes(con) == 2 ? true : false
is_hyperlinkconstr(con::LinkConstraint) = getnumnodes(con) > 2 ? true : false

#Extend JuMP's add constraint for link models.  Return a reference to the constraint
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
