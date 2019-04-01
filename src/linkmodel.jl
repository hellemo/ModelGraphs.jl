#####################################################
# Link Model
# A link model stores model data about graph links
#####################################################
mutable struct LinkModel <: JuMP.AbstractModel   #subtyping here so I can get ConstraintRef

    graphvariables::Dict{Int,GraphVariable}                             #global level variables.  Defined over an entire graph.  Can be used by sub-problems.
    graphvariablechildren::Dict{GraphVariable,JuMP.AbstractVariable}    #children variables of a graph variable
    graphconstraints::Dict{Int,JuMP.AbstractConstraint}  #global constraint.  Defined over graph variables.  Think "first stage" constraints or Shared variable constraints.
    linkconstraints::Dict{Int,LinkConstraint}            #linking constraint.  Defined over variables in nodes.
    allconstraints::Dict{Int,JuMP.AbstractConstraint}    #Both graph and link constraints

    graphvarnames::Dict{Int,String}
    graphconstraintnames::Dict{Int,String}
    linkconstraintnames::Dict{Int,String}

    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar
    solver::Union{MOI.AbstractOptimizer,AbstractGraphSolver}

    graph_var_index::Int          #Track indices so we can manage deleting variables and constraints
    graph_constraint_index::Int
end

#Constructor
#LinkModel(;solver = MOI.) = LinkModel(LinkData(),0,JuMP.AffExpr(),solver)
# model = new(
#                    0, Dict{Int, JuMP.AbstractVariable}(),   Dict{Int, String}(),    # Model Variables
#                    0, Dict{Int, JuMP.AbstractConstraint}(), Dict{Int, String}(),    # Model Constraints
#                    MOI.FEASIBILITY_SENSE, zero(JuMP.GenericAffExpr{Float64, StructuredVariableRef}), # Model objective
# Dict{Symbol, Any}())

# Graph Constraint Reference
struct GraphConstraintRef
    model::LinkModel # `model` owning the constraint
    idx::Int         #  index in `model.graphconstraints`
end
JuMP.constraint_type(::LinkModel) = GraphConstraintRef
#####################################################
#  Link Constraint
#  A linear constraint between JuMP Models (nodes)
#  Link constraints can be equality or inequality.
#####################################################
#abstract type LinkSet <: MOI.AbstractScalarSet end
struct LinkConstraint{F <: AbstractJuMPScalar,S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::F
    set::S
    node_indices::Vector{Int}  #Custom link constraint data.  Indices of nodes this LinkConstraint connects.

    function LinkConstraint(func::F,set::S)
        linkcon = new()
        linkcon.func = func
        linkcon.set = set
        linkcon.node_indices = [getnode(var) for var in keys(func.terms)]
    end
end

function addlinkconstraint(m::LinkModel, con::LinkConstraint, name::String="")
    m.link_constraint_index += 1
    cref = GraphConstraintRef(m, m.link_constraint_index)
    m.linkconstraints[cref.idx] = linkcon
    JuMP.set_name(cref, name)
    return cref
end

# function JuMP.add_constraint(m::LinkModel, con::JuMP.AbstractConstraint, name::String="")
#     error("ModelGraphs only support Scalar Link Constraints.")
# end

#Add a LinkConstraint to a LinkModel
# function JuMP.add_constraint(m::LinkModel, con::LinkConstraint, name::String="")
#     m.link_constraint_index += 1
#     cref = GraphConstraintRef(m, m.link_constraint_index)
#     m.linkconstraints[cref.idx] = linkcon
#     JuMP.set_name(cref, name)
#     return cref
# end

#NOTE build_constraint will create a standard linear constraint.  I will then convert this to a LinkConstraint in the @linkconstraint macro
#Build a Link Constraint.  A Link Constraint is an AffineExpression defined over a LinkSet (GreatThan, LessThan, EqualTo)
# function JuMP.build_constraint(_error::Function, func::JuMP.AffExpr,set::LinkSet)#,nodes::Vector{Int})
#     nodes = [getnode(var) for var in keys(func.terms)]
#     link_constraint = LinkConstraint(func,set,nodes)
#     return link_constraint
# end
jump_function(constraint::LinkConstraint) = constraint.func
moi_set(constraint::LinkConstraint) = constraint.set
shape(::LinkConstraint) = JuMP.ScalarShape()

#TODO Figuring out how this will work
# function constraint_object(ref::ConstraintRef{Model, _MOICON{FuncType, LinkSet}}) where {FuncType <: MOI.AbstractScalarFunction}
#     model = ref.model
#     f = MOI.get(model, MOI.ConstraintFunction(), ref)::FuncType
#     s = MOI.get(model, MOI.ConstraintSet(), ref)::LinkSet
#     return ScalarConstraint(jump_function(model, f), s)
# end
function constraint_object(cref::GraphConstraintRef, F::Type, S::Type}
   con = cref.model.allconstraints[cref.idx]
   con.func::F
   con.set::S
   return con
end

function MOI.delete!(m::LinkModel, cref::GraphConstraintRef)
    if typeof(constraint_object(cref)) == LinkConstraint
        delete!(m.linkconstraints, cref.idx)
        delete!(m.linkconnames, cref.idx)
    elseif typeof(constraint_object(cref)) == JuMP.AbstractConstraint
        delete!(m.graphconstraints, cref.idx)
        delete!(m.graphconnames, cref.idx)
    else
        error("Constraint reference refers to constraint type that was not recognized")
    end
end

MOI.is_valid(m::LinkModel, cref::GraphConstraintRef) = cref.idx in keys(m.allconstraints)
# function JuMP.constraint_object(cref::GraphConstraintRef, F::Type, S::Type)
#     con = cref.model.graphconstraints[cref.idx]
#     con.func::F
#     con.set::S
#     return con
# end

# #I think this will work?  Might need to do if S == LinkSet or something
# function JuMP.constraint_object(cref::GraphConstraintRef, F::Type, S::LinkSet)
#     con = cref.model.linkconstraints[cref.idx]
#     con.func::F
#     con.set::S
#     return con
# end

#####################################################
#  Graph Variables
#  A graph variable belongs to the graph as opposed to an individual model.  A graph variable is available for
#  different models to use in their constraints.  This implies it is a shared variable between nodes in the graph.
#####################################################
#Need this to dispatch @constraint with graph variables?
# NOTE: Not sure I need this if I have the GraphVariableRef to dispatch on
# struct GraphVariable <: JuMP.AbstractVariable
#     info::JuMP.VariableInfo{S, T, U, V}
# end

#Reference to a Graph Variable.  These variables can be referenced in constraints from JuMP models
struct GraphVariableRef <: JuMP.AbstractVariableRef  #NOTE AbstractVariableRef is an AbstractJuMPScalar
    model::LinkModel
    idx::Int
end

Base.broadcastable(v::GraphVariableRef) = Ref(v)
Base.copy(v::GraphVariableRef) = v
Base.:(==)(v::GraphVariableRef, w::GraphVariableRef) = v.model === w.model && v.idx == w.idx
JuMP.owner_model(v::GraphVariableRef) = v.model
JuMP.isequal_canonical(v::GraphVariableRef, w::GraphVariableRef) = v == w
JuMP.variable_type(::LinkModel) = GraphVariableRef

function JuMP.add_variable(m::LinkModel, v::JuMP.AbstractVariable, name::String="")
    m.graph_var_index += 1
    vref = GraphVariableRef(m, m.graph_var_index)
    m.graphvariables[vref.idx] = v
    JuMP.set_name(vref, name)
    return vref
end

function MOI.delete!(m::LinkModel, vref::GraphVariableRef)
    delete!(m.linkvariables, vref.idx)
    delete!(m.linkvarnames, vref.idx)
end
MOI.is_valid(m::LinkModel, vref::GraphVariableRef) = vref.idx in keys(m.linkvariables)
JuMP.num_variables(m::LinkModel) = length(m.linkvariables)

#####################################################
#  Graph Constraint
#  A constraint between graph variables. Think first-stage constraints.
#####################################################
# const GraphAffExpr = JuMP.GenericAffExpr{AffExpr,GraphVariableRef}
# function JuMP.build_constraint(_error::Function,func::GraphAffExpr,set::MOI.AbstractScalarSet)
#     constraint = ScalarConstraint(func, set)
#     return constraint
# end

#Add a GraphConstraint to a LinkModel
function JuMP.add_constraint(m::LinkModel, con::AbstractConstraint, name::String="")
    m.graph_constraint_index += 1
    cref = GraphConstraintRef(m, m.graph_constraint_index)
    m.graphconstraints[cref.idx] = con
    JuMP.set_name(cref, name)
    return cref
end


# Add constraints on nodes that use a GraphVariableRef
function JuMP.add_constraint(m::Model, con::AbstractConstraint, name::String="")
    m.graph_constraint_index += 1
    cref = GraphConstraintRef(m, m.graph_constraint_index)
    m.graphconstraints[cref.idx] = con
    JuMP.set_name(cref, name)
    return cref
end

#####################################################
#  Shared Variable Constraint
#  A constraint on a model that uses Graph Variables
#####################################################
#func::GraphVariableRef
#TODO: I think I need a custom variable type to dispatch on @constraint for node models.  Node models could add constraints that use the graph model.



########################################################################
# Node information for graph and link constraints
########################################################################
function StructureGraphs.getnodes(con::LinkConstraint)
    # vars = keys(con.func.terms)
    # nodes = unique([getnode(var) for var in vars])  #NOTE unique would lead to slow behavior
    # return nodes
    return [getnode(con.nodes[i]) for i = 1:getnumnodes(con)]
end

getnumnodes(con::LinkConstraint) = return length(con.nodes)



#getlinkdata(model::LinkModel) = model.linkdata

#Get the 2 variable or multi-variable linkconstraints
#getlinkconstraints(model::LinkModel) = getlinkdata(model).linkconstraints
#getsimplelinkconstraints(model::LinkModel) = getlinkdata(model).linkconstraints[getlinkdata(model).simple_links]
#gethyperlinkconstraints(model::LinkModel) = getlinkdata(model).linkconstraints[getlinkdata(model).hyper_links]

is_simplelinkconstr(con::LinkConstraint) = getnumnodes(con) == 2 ? true : false
is_hyperlinkconstr(con::LinkConstraint) = getnumnodes(con) > 2 ? true : false

#Extend JuMP's add constraint for link models.  Return a reference to the constraint
# #NOTE JuMP no longer uses LinearConstraint
# function JuMP.addconstraint(model::LinkModel,constr::JuMP.LinearConstraint)
#     #TODO Do some error checking here
#     linkdata = getlinkdata(model)
#     linkconstr = LinkConstraint(constr)
#     push!(linkdata.linkconstraints,linkconstr)
#     ref = ConstraintRef{LinkModel,LinkConstraint}(model, length(linkdata.linkconstraints))
#
#     if is_linkconstr(linkconstr)
#         push!(linkdata.simple_links,length(linkdata.linkconstraints))
#     elseif is_hyperconstr(linkconstr )
#         push!(linkdata.hyper_links,length(linkdata.linkconstraints))
#     else
#         error("constraint $constr doesn't make sense")
#     end
#     return ref
# end


#Constructor
#NOTE Likely not used anymore
# LinkConstraint(con::JuMP.LinearConstraint) = LinkConstraint(con.terms,con.lb,con.ub)
# LinkConstraint(ref::ConstraintRef) = ref.m.linkdata.linkconstraints[ref.idx]  #Get the Link constraint from a constraint reference
# LinkConstraint(con::LinkConstraint) = LinkConstraint(con.terms,con.lb,con.ub)
