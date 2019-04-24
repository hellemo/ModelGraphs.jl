#####################################################
# Link Model
# A link model stores model data about graph links
#####################################################
mutable struct LinkModel <:  AbstractLinkModel
    graph::AbstractModelGraph

    graphvariables::Dict{Int,JuMP.AbstractVariable}                        #global level variables.  Defined over an entire graph.  Can be used by sub-problems.
    graphvariablemap::Dict{JuMP.AbstractVariable,JuMP.AbstractVariable}    #map of graph variables to children variables

    graphconstraints::Dict{Int,JuMP.AbstractConstraint}                    #global constraint.  Defined over graph variables.  Think "first stage" constraints or Shared variable constraints.
    linkconstraints::Dict{Int,AbstractLinkConstraint}                      #linking constraint.  Defined over variables in nodes.

    graphvarnames::Dict{Int,String}
    graphconstraintnames::Dict{Int,String}
    linkconstraintnames::Dict{Int,String}

    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar
    #solver::Union{MOI.AbstractOptimizer,AbstractGraphSolver,Nothing}  #NOTE: Might use a JuMP SolverFactory here

    graph_var_index::Int                                                    #Track indices so we can manage deleting variables and constraints
    graph_constraint_index::Int

    objdict::Dict{Symbol,Any}

    #TODO nlp_data.  GraphConstraints could be nonlinear.

    #Constructor
    function LinkModel(graph::AbstractModelGraph)
        model = new(graph,
                    Dict{Int, JuMP.AbstractVariable}(),
                    Dict{JuMP.AbstractVariable, JuMP.AbstractVariable}(),
                    Dict{Int, JuMP.AbstractConstraint}(),
                    Dict{Int, AbstractLinkConstraint}(),
                    Dict{Int, String}(),
                    Dict{Int, String}(),
                    Dict{Int, String}(),
                    MOI.FEASIBILITY_SENSE,
                    zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef}),
                    #nothing,
                    0,
                    0,
                    Dict{Symbol,Any}()
                    )
    end
end
JuMP.object_dictionary(m::LinkModel) = m.objdict
JuMP.objective_sense(m::LinkModel) = m.objective_sense
getgraphvariables(m::LinkModel) = collect(values(m.graphvariables))
getgraphconstraints(m::LinkModel) = collect(values(m.graphconstraints))
getlinkconstraints(model::LinkModel) = collect(values(model.linkconstraints))
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
struct GraphConstraintRef <: AbstractGraphConstraintRef
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

#####################################################
#  Link Constraint
#  A linear constraint between JuMP Models (nodes).  Link constraints can be equality or inequality.
#####################################################
struct LinkConstraint{F <: JuMP.AbstractJuMPScalar,S <: MOI.AbstractScalarSet} <: AbstractLinkConstraint
    func::F
    set::S
    #NOTE:
    #Could cache other data on the LinkConstraint
end

LinkConstraint(con::JuMP.ScalarConstraint) = LinkConstraint(con.func,con.set)
LinkConstraint(ref::GraphConstraintRef) = JuMP.owner_model(ref).linkconstraints[ref.idx]

#Add a LinkConstraint to a LinkModel
function JuMP.add_constraint(m::LinkModel, con::JuMP.ScalarConstraint, name::String="")
    m.graph_constraint_index += 1
    cref = GraphConstraintRef(m, m.graph_constraint_index)
    link_con = LinkConstraint(con)      #convert ScalarConstraint to a LinkConstraint
    m.linkconstraints[cref.idx] = link_con
    JuMP.set_name(cref, name)

    #Add LinkingEdges to the Graph
    graph = m.graph
    addlinkedges!(graph,cref)

    return cref
end

function JuMP.add_constraint(m::LinkModel, con::JuMP.AbstractConstraint, name::String="")
    error("Link Models only support graph constraints and link constraints")
end

jump_function(constraint::LinkConstraint) = constraint.func
moi_set(constraint::LinkConstraint) = constraint.set
shape(::LinkConstraint) = JuMP.ScalarShape()
