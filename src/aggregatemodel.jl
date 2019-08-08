"""
    AggregateGraph
    Extension Data attached to an Aggregated ModelGraph.  The AggregateGraph retains a reference to the original ModelGraph topology with references to
    graphvariables, graphconstraints, and linkconstraints.
"""
#Aggregation Data
#After aggregating, we throw away the topology, but we keep references to what was a linkvariable and linkconstraint
mutable struct AggregationInfo # <: AbstractModelGraph
    #hypergraph::NestedHyperGraph
    nodes::Vector{NodeReference}
    linkvariables::Vector{VariableRef}
    linkconstraints::Vector{ConstraintRef}
    NLlinkconstraint::Vector{ConstraintRef}
end
AggregateInfo() = AggregateModel()
gethypergraph(graph::AggregateGraph) = graph.hypergraph

mutable struct NodeReference #<: AbstractModelNode
    hypernode::HyperNode
    obj_dict::Dict{Symbol,Any}
    variablemap::Dict{JuMP.VariableRef,JuMP.VariableRef}
    constraintmap::Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}
    nl_constraintmap::Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}
    objective::Union{JuMP.AbstractJuMPScalar,Expr}
end

zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef}))

#Construct a structured model, but roll it all into one JuMP model (this is how we solve with JuMP accessible solvers)
function AggregateModel()
    m = JuMP.Model()
    m.ext[:AggregationInfo] = AggregationInfo()
    return m
end

is_aggregate_model(m::JuMP.Model) = haskey(m.ext,:AggregationInfo) ? true : false  #check if the model is a graph model
assert_aggregate_model(m::JuMP.Model) = @assert is_aggregate_model(m)
#Define all of the JuMP model extension functions
get_aggregation_info(m::JuMP.Model) = haskey(m.ext, :AggregationInfo) ? m.ext[:AggregationInfo] : error("Model is not a graph model")


JuMP.objective_function(node::NodeReference) = node.objective

getnodevariables(node::NodeReference) = node.variablelist
getnodevariable(node::NodeReference,index::Integer) = node.variablelist[index]
getnodeconstraints(node::NodeReference) = node.constraintlist
JuMP.num_variables(node::NodeReference) = length(node.variablelist)
#get node variables using a symbol lookup
getindex(node::NodeReference,s::Symbol) = node.obj_dict[s]

#get all of the link constraints from a JuMP model
#Retrieves constraint indices that are link constraints
function get_link_constraints(m::JuMP.Model)
    assert_aggregate_model(m)
    agg_info = get_aggregation_info(m)
    return agg_info.linkconstraints
end


# Should be defined by base type
# #Add nodes and edges to graph models.  These are used for model instantiation from a graph
# function NestedHyperGraphs.add_node!(m::JuMP.Model; index = nv(getgraph(m).hypergraph)+1)
#     is_graphmodel(m) || error("Can only add nodes to graph models")
#     node = create_node(getgraph(m))
#     add_node!(getgraph(m),node,index = index)
#     return node
# end
#
# function StructureGraphs.add_edge!(m::JuMP.Model,nodes::NodeReference...)
#     is_graphmodel(m) || error("Can only add edges to graph models")
#     edge = add_edge!(getgraph(m),nodes...)
#     return edge
# end

# StructureGraphs.getnodes(m::JuMP.Model) = getnodes(getgraph(m))
# StructureGraphs.getedges(m::JuMP.Model) = getedges(getgraph(m))
#
# StructureGraphs.getnode(m::JuMP.Model,id::Integer) = getnode(getgraph(m),id)  #Grab from the highest level graph if not specified
# StructureGraphs.getnode(m::JuMP.Model,sid::Integer,nid::Integer) = getnode(getgraph(m).subgraphlist[sid],nid)

# StructureGraphs.getedge(m::Model,id::LightGraphs.AbstractEdge) = getedge(getgraph(m))[id]
# StructureGraphs.getedge(m::Model,sid::Integer,eid::LightGraphs.AbstractEdge) = getedge(getgraph(m).subgraphlist[sid],eid)
