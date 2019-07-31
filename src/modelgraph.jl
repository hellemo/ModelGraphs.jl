##############################################################################
# ModelGraph
##############################################################################
"""
ModelGraph()

The ModelGraph Type.  Represents a graph containing models (nodes) and the linkconstraints (edges) between them.
A ModelGraph wraps a BasePlasmoGraph and can use its methods.  A ModelGraph also wraps a LinkModel object which extends a JuMP AbstractModel to provide model management functions.

"""
mutable struct ModelGraph <: AbstractModelGraph
    hypergraph::NestedHyperGraph                                              #This is a Multilevel Hypergraph

    mastermodel::JuMP.Model                                               #JuMP Model we use to store link variables and master constraints
    linkvariables::Dict{Int,JuMP.AbstractVariable}                       #Link Variable reference.  These variables are also in the mastermodel
    linkvariablemap::Dict{JuMP.AbstractVariable,JuMP.AbstractVariable}    #map of link variables in master model to corresponding variables in ModelNodes
    linkvariablenames::Dict{Int,String}


    linkconstraints::Dict{Int,AbstractLinkConstraint}                      #linking constraint.  Defined over variables in ModelNodes.
    linkconstraintnames::Dict{Int,String}
    linkconstraint_hyperedge_map::Dict{AbstractLinkConstraint,HyperEdge}

    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar
    optimizer::Union{JuMP.OptimizerFactory,AbstractGraphSolver,Nothing}

    link_variable_index::Int                                                    #Track indices so we can manage deleting variables and constraints
    graph_constraint_index::Int  #keep track of master and link constraints

    # link_constraint_index::Int
    # master_constraint_index::Int


    objdict::Dict{Symbol,Any}

    nlp_data::JuMP._NLPData                                                 #This is strictly for nonlinear-link constraints

    #TODO nlp_data? : Using the mastermodel for NLP data

    #Constructor
    function ModelGraph()
        model = new(
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
        return model
    end
end

#HyperGraph Interface
gethypergraph(graph::ModelGraph) = graph.hypergraph

function add_subgraph!(graph::ModelGraph,subgraph::ModelGraph)
    hypergraph = gethypergraph(graph)
    sub_hypergraph = gethypergraph(subgraph)
    add_subgraph!(hypergraph,subgraph)
    return true
end



getnumlinkconstraints(graph::AbstractModelGraph) = length(graph.linkconstraints)
getnumnllinkconstraint(graph::AbstractModelGraph) = "not supported yet"

hasobjective(graph::AbstractModelGraph) = getlinkmodel(graph).objective_function != zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef})

"Set the objective of a ModelGraph"
#TODO. Write more objective methods for the LinkModel
JuMP.set_objective_function(graph::AbstractModelGraph, sense::MOI.OptimizationSense, x::JuMP.VariableRef) = JuMP.set_objective_function(graph, sense, convert(AffExpr,x))

"Get the ModelGraph objective value"
JuMP.getobjectivevalue(graph::AbstractModelGraph) = getobjectivevalue(graph.linkmodel)

"Get the current created JuMP model for the ModelGraph.  Only created when solving using a JuMP compliant solver."
getinternaljumpmodel(graph::AbstractModelGraph) = graph.serial_model

getlinkvariables(graph::AbstractModelGraph) = graph.link
getmasterconstraints(graph::AbstractModelGraph) = getgraphconstraints(getlinkmodel(graph))

"""
get_all_linkconstraints(graph::AbstractModelGraph)

Get a list containing every link constraint in the graph, including its subgraphs
"""
function get_all_linkconstraints(graph::AbstractModelGraph)
    links = []
    for subgraph in getsubgraphlist(graph)
        append!(links,getlinkconstraints(subgraph))
    end
    append!(links,getlinkconstraints(graph))
    return links
end

JuMP.object_dictionary(m::ModelGraph) = m.objdict
JuMP.objective_sense(m::ModelGraph) = m.objective_sense
getlinkvariables(m::ModelGraph) = collect(values(m.graphvariables))
getmasterconstraints(m::ModelGraph) = collect(values(m.graphconstraints))
getlinkconstraints(model::ModelGraph) = collect(values(model.linkconstraints))
#####################################################
#  Link Variables
#  A graph variable belongs to the graph as opposed to an individual model.  A graph variable is available for
#  different models to use in their constraints.  This implies it is a shared variable between nodes in the graph.
#####################################################
#Reference to a Link Variable.  These variables can be referenced in constraints from JuMP models
struct LinkVariableRef <: JuMP.AbstractVariableRef  #NOTE AbstractVariableRef is an AbstractJuMPScalar
    model::ModelGraph
    idx::Int                    #index in the master model
end

Base.broadcastable(v::LinkVariableRef) = Ref(v)
Base.copy(v::LinkVariableRef) = v
Base.:(==)(v::LinkVariableRef, w::LinkVariableRef) = v.model === w.model && v.idx == w.idx
JuMP.owner_model(v::LinkVariableRef) = v.model
JuMP.isequal_canonical(v::LinkVariableRef, w::LinkVariableRef) = v == w
JuMP.variable_type(::LinkModel) = LinkVariableRef

JuMP.set_name(v::LinkVariableRef, s::String) = JuMP.owner_model(v).graphvarnames[v.idx] = s
JuMP.name(v::LinkVariableRef) =  JuMP.owner_model(v).graphvarnames[v.idx]

function JuMP.add_variable(m::ModelGraph, v::JuMP.AbstractVariable, name::String="")
    m.link_var_index += 1
    vref = LinkVariableRef(m, m.link_var_index)
    JuMP.add_variable(m.mastermodel,v,name)  #add to master model
    m.linkvariables[vref.idx] = v
    JuMP.set_name(vref, name)
    return vref
end

function MOI.delete!(m::ModelGraph, vref::LinkVariableRef)
    delete!(m.graphvariables, vref.idx)
    delete!(m.graphvarnames, vref.idx)
end
MOI.is_valid(m::LinkModel, vref::LinkVariableRef) = vref.idx in keys(m.linkvariables)
JuMP.num_variables(m::ModelGraph) = length(m.linkvariables)

#####################################################
# Graph Constraint (Master Constraint)
# A constraint between graph variables. Think first-stage constraints.
#####################################################
const GraphAffExpr = JuMP.GenericAffExpr{Float64,LinkVariableRef}

# Constraint Reference to a Master or LinkConstraint
struct GraphConstraintRef <: AbstractGraphConstraintRef
    modelgraph::ModelGraph #       `model` owning the constraint
    idx::Int                       #index in `model.graphconstraints`
end
JuMP.constraint_type(::ModelGraph) = GraphConstraintRef
JuMP.owner_model(con::GraphConstraintRef) = con.modelgraph

function JuMP.constraint_object(cref::GraphConstraintRef, F::Type, S::Type)
   con = cref.modelgraph.allconstraints[cref.idx]
   con.func::F
   con.set::S
   return con
end

function MOI.delete!(m::ModelGraph, cref::GraphConstraintRef)
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

MOI.is_valid(m::ModelGraph, cref::GraphConstraintRef) = cref.idx in keys(m.)#keys(m.graphconstraints) || cref.idx in keys(m.linkconstraints)

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

#Add a Master Constraint
function JuMP.add_constraint(m::ModelGraph, con::ScalarGraphConstraint, name::String="")
    m.graph_constraint_index += 1
    cref = GraphConstraintRef(m, m.graph_constraint_index)
    JuMP.add_constraint(m.mastermodel,con,name)   #also add to master model
    m.graphconstraints[cref.idx] = con
    JuMP.set_name(cref, name)
    return cref
end


# Model Extras
JuMP.show_constraints_summary(::IOContext,m::ModelGraph) = ""
JuMP.show_backend_summary(::IOContext,m::ModelGraph) = ""

#####################################################
#  Link Constraint
#  A linear constraint between JuMP Models (nodes).  Link constraints can be equality or inequality.
#####################################################
struct LinkConstraint{F <: JuMP.AbstractJuMPScalar,S <: MOI.AbstractScalarSet} <: AbstractLinkConstraint
    func::F
    set::S
    # #Caching data on the LinkConstraint
    # graph::AbstractModelGraph
    #node_indices::Vector{Int64}   #Should be easier to just reference a hyperedge
    #hyperedge::HyperEdge  #reference to a hyperedge
end
LinkConstraint(ref::GraphConstraintRef) = JuMP.owner_model(ref).linkconstraints[ref.idx]

# function LinkConstraint(con::JuMP.ScalarConstraint,graph::AbstractModelGraph)
#     node_indices = sort(unique([getindex(graph,getnode(var)) for var in keys(con.func.terms)]))
#     return LinkConstraint(con.func,con.set,graph,node_indices)
# end
function LinkConstraint(con::JuMP.ScalarConstraint,hyperedge::GraphEdge)
    #node_indices = sort(unique([getindex(graph,getnode(var)) for var in keys(con.func.terms)]))
    return LinkConstraint(con.func,con.set,hyperedge)
end

function getnodes(con::LinkConstraint)
    #TODO: Check uniqueness.  It should be unique now that JuMP uses an OrderedDict to store terms.
    #return [getnode(var) for var in keys(con.func.terms)]
    return getnodes(con.hyperedge)
    #return map(n -> getnode(con.graph,n),con.node_indices)
end
getnumnodes(con::LinkConstraint) = length(getnodes(con))

#Add a LinkConstraint to a LinkModel
function JuMP.add_constraint(m::ModelGraph, con::JuMP.ScalarConstraint, name::String="")
    m.graph_constraint_index += 1
    cref = GraphConstraintRef(m, m.graph_constraint_index)

    #Setup graph information
    hypergraph = gethypergraph(m)
    node_indices = sort(unique([getindex(hypergraph,getnode(var)) for var in keys(con.func.terms)]))
    hyperedge = getedge(hypergraph,node_indices)

    link_con = LinkConstraint(con,graph,hyperedge)      #convert ScalarConstraint to a LinkConstraint
    m.linkconstraints[cref.idx] = link_con
    JuMP.set_name(cref, name)

    #Add LinkingEdges to the Graph
    #addlinkedges!(graph,cref)
    add_linkconstraint_edges!(hypergraph,hyperedge)

    return cref
end

function JuMP.add_constraint(m::LinkModel, con::JuMP.AbstractConstraint, name::String="")
    error("A ModelGraph only supports LinkConstraints")
end

jump_function(constraint::LinkConstraint) = constraint.func
moi_set(constraint::LinkConstraint) = constraint.set
shape(::LinkConstraint) = JuMP.ScalarShape()


#################################
# Solver setters and getters
#################################
#TODO Figure out how this works in JuMP 0.19.  Might have to hook into MOI here.
#set_optimizer(model::AbstractModelGraph,graph_solver::AbstractGraphSolver) = set_optimizer(model.linkmodel,graph_solver)

####################################
#Print Functions
####################################
function string(graph::ModelGraph)
    "Model Graph:"*string(getlabel(graph))*"
nodes:"*string((length(getnodes(graph))))*"
link constraints (edges):"*string((getnumlinkconstraints(graph)))
end
print(io::IO, graph::AbstractModelGraph) = print(io, string(graph))
show(io::IO,graph::AbstractModelGraph) = print(io,graph)
