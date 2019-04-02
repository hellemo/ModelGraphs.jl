#abstract type AbstractGraphVariableRef <: JuMP.AbstractVariableRef
abstract type AbstractLinkConstraint <: JuMP.AbstractConstraint end
#####################################################
# Link Model
# A link model stores model data about graph links
#####################################################
mutable struct LinkModel <: JuMP.AbstractModel   #subtyping here so I can get ConstraintRef

    graphvariables::Dict{Int,JuMP.AbstractVariable}                             #global level variables.  Defined over an entire graph.  Can be used by sub-problems.
    graphvariablechildren::Dict{JuMP.AbstractVariable,JuMP.AbstractVariable}    #map of graph variables to children variables

    graphconstraints::Dict{Int,JuMP.AbstractConstraint}  #global constraint.  Defined over graph variables.  Think "first stage" constraints or Shared variable constraints.
    linkconstraints::Dict{Int,AbstractLinkConstraint}            #linking constraint.  Defined over variables in nodes.

    graphvarnames::Dict{Int,String}
    graphconstraintnames::Dict{Int,String}
    linkconstraintnames::Dict{Int,String}

    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar
    solver::Union{MOI.AbstractOptimizer,AbstractGraphSolver,Nothing}

    graph_var_index::Int          #Track indices so we can manage deleting variables and constraints
    graph_constraint_index::Int

    objdict::Dict{Symbol,Any}

    #Constructor
    function LinkModel()
        model = new(Dict{Int, JuMP.AbstractVariable}(),
                    Dict{JuMP.AbstractVariable, JuMP.AbstractVariable}(),
                    Dict{Int, JuMP.AbstractConstraint}(),
                    Dict{Int, AbstractLinkConstraint}(),
                    Dict{Int, String}(),
                    Dict{Int, String}(),
                    Dict{Int, String}(),
                    MOI.FEASIBILITY_SENSE,
                    zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef}),
                    nothing,
                    0,
                    0,
                    Dict{Symbol,Any}()
                    )
    end
end
JuMP.object_dictionary(m::LinkModel) = m.objdict
JuMP.objective_sense(m::LinkModel) = m.objective_sense

#####################################################
#  Graph Variables
#  A graph variable belongs to the graph as opposed to an individual model.  A graph variable is available for
#  different models to use in their constraints.  This implies it is a shared variable between nodes in the graph.
#####################################################
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

JuMP.set_name(v::GraphVariableRef, s::String) = JuMP.owner_model(v).graphvarnames[v.idx] = s
JuMP.name(v::GraphVariableRef) =  JuMP.owner_model(v).graphvarnames[v.idx]

function JuMP.add_variable(m::LinkModel, v::JuMP.AbstractVariable, name::String="")
    m.graph_var_index += 1
    vref = GraphVariableRef(m, m.graph_var_index)
    m.graphvariables[vref.idx] = v
    JuMP.set_name(vref, name)
    return vref
end

function MOI.delete!(m::LinkModel, vref::GraphVariableRef)
    delete!(m.graphvariables, vref.idx)
    delete!(m.graphvarnames, vref.idx)
end
MOI.is_valid(m::LinkModel, vref::GraphVariableRef) = vref.idx in keys(m.graphvariables)
JuMP.num_variables(m::LinkModel) = length(m.graphvariables)

#####################################################
# Graph Constraint
# A constraint between graph variables. Think first-stage constraints.
# TODO Extend @constraint macro to work with graph variables
#####################################################
const GraphAffExpr = JuMP.GenericAffExpr{Float64,GraphVariableRef}

# Graph Constraint Reference
struct GraphConstraintRef
    model::LinkModel # `model` owning the constraint
    idx::Int         #  index in `model.graphconstraints`
end
JuMP.constraint_type(::LinkModel) = GraphConstraintRef
JuMP.owner_model(con::GraphConstraintRef) = con.model

function constraint_object(cref::GraphConstraintRef, F::Type, S::Type)
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

MOI.is_valid(m::LinkModel, cref::GraphConstraintRef) = cref.idx in keys(m.graphconstraints) || cref.idx in keys(m.linkconstraints)

JuMP.set_name(con::GraphConstraintRef, s::String) = JuMP.owner_model(con).graphconstraintnames[con.idx] = s
JuMP.name(con::GraphConstraintRef) =  JuMP.owner_model(con).graphconstraintnames[con.idx]

struct ScalarGraphConstraint{F <: GraphAffExpr,S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::F
    set::S
    function ScalarGraphConstraint(F::AbstractJuMPScalar,S::MOI.AbstractScalarSet)
        con = new{typeof(F),typeof(S)}(F,S)
    end
end

function JuMP.build_constraint(_error::Function,func::GraphAffExpr,set::MOI.AbstractScalarSet)
    constraint = ScalarGraphConstraint(func, set)
    return constraint
end

#Fallback.  Add any Constraint to a LinkModel
function JuMP.add_constraint(m::LinkModel, con::ScalarGraphConstraint, name::String="")
    m.graph_constraint_index += 1
    cref = GraphConstraintRef(m, m.graph_constraint_index)
    m.graphconstraints[cref.idx] = con
    JuMP.set_name(cref, name)
    return cref
end

# Model Extras
JuMP.show_constraints_summary(::IOContext,m::LinkModel) = ""
JuMP.show_backend_summary(::IOContext,m::LinkModel) = ""
########################################################################
# Node information for graph and link constraints
########################################################################
# function StructureGraphs.getnodes(con::LinkConstraint)
#     # vars = keys(con.func.terms)
#     # nodes = unique([getnode(var) for var in vars])  #NOTE unique would lead to slow behavior
#     # return nodes
#     return [getnode(con.nodes[i]) for i = 1:getnumnodes(con)]
# end
# getnumnodes(con::LinkConstraint) = return length(con.nodes)

# Add a graph constraint to a standard JuMP model
# # TODO
# function JuMP.add_constraint(m::Model, con::GraphConstraint, name::String="")
#     # m.graph_constraint_index += 1
#     # cref = GraphConstraintRef(m, m.graph_constraint_index)
#     # m.graphconstraints[cref.idx] = con
#     # JuMP.set_name(cref, name)
#     # return cref
# end




# function addlinkconstraint(m::LinkModel, linkcon::LinkConstraint, name::String="")
#     m.graph_constraint_index += 1
#     cref = GraphConstraintRef(m, m.link_constraint_index)
#     m.linkconstraints[cref.idx] = linkcon
#     JuMP.set_name(cref, name)
#     return cref
# end







# abstract type LinkSet <: MOI.AbstractScalarSet end




#Need this to dispatch @constraint with graph variables?
# NOTE: Not sure I need this if I have the GraphVariableRef to dispatch on
# struct GraphVariable <: JuMP.AbstractVariable
#     info::JuMP.VariableInfo{S, T, U, V}
# end




#####################################################
#  Shared Variable Constraint
#  A constraint on a model that uses Graph Variables
#####################################################
#func::GraphVariableRef
#TODO: I think I need a custom variable type to dispatch on @constraint for node models.  Node models could add constraints that use the graph model.







#getlinkdata(model::LinkModel) = model.linkdata

#Get the 2 variable or multi-variable linkconstraints
#getlinkconstraints(model::LinkModel) = getlinkdata(model).linkconstraints
#getsimplelinkconstraints(model::LinkModel) = getlinkdata(model).linkconstraints[getlinkdata(model).simple_links]
#gethyperlinkconstraints(model::LinkModel) = getlinkdata(model).linkconstraints[getlinkdata(model).hyper_links]



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

#NOTE build_constraint will create a standard linear constraint.  I will then convert this to a LinkConstraint in the @linkconstraint macro
#Build a Link Constraint.  A Link Constraint is an AffineExpression defined over a LinkSet (GreatThan, LessThan, EqualTo)
# function JuMP.build_constraint(_error::Function, func::JuMP.AffExpr,set::LinkSet)#,nodes::Vector{Int})
#     nodes = [getnode(var) for var in keys(func.terms)]
#     link_constraint = LinkConstraint(func,set,nodes)
#     return link_constraint
# end

#Constructor
#NOTE Likely not used anymore
# LinkConstraint(con::JuMP.LinearConstraint) = LinkConstraint(con.terms,con.lb,con.ub)
# LinkConstraint(ref::ConstraintRef) = ref.m.linkdata.linkconstraints[ref.idx]  #Get the Link constraint from a constraint reference
# LinkConstraint(con::LinkConstraint) = LinkConstraint(con.terms,con.lb,con.ub)

#TODO Figuring out how this will work
# function constraint_object(ref::ConstraintRef{Model, _MOICON{FuncType, LinkSet}}) where {FuncType <: MOI.AbstractScalarFunction}
#     model = ref.model
#     f = MOI.get(model, MOI.ConstraintFunction(), ref)::FuncType
#     s = MOI.get(model, MOI.ConstraintSet(), ref)::LinkSet
#     return ScalarConstraint(jump_function(model, f), s)
# end
