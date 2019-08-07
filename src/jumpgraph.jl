"""
    AggregateGraph
    Extension Data attached to an Aggregated ModelGraph.  The AggregateGraph retains a reference to the original ModelGraph topology with references to
    graphvariables, graphconstraints, and linkconstraints.
"""
mutable struct AggregateGraph <: AbstractModelGraph
    hypergraph::NestedHyperGraph
    master::JuMP.Model
    linkvariables::Vector{LinkVariableRef}
    linkconstraints::Vector{LinkConstraintRef}
end
AggregateModel() = AggregateModel(StructureGraphs.StructuredHyperGraph(),GraphVariableRef[],GraphConstraintRef[],GraphConstraintRef[])
JuMPGraph(hypergraph::StructureGraphs.StructureGraph) = JuMPGraph(hypergraph,GraphVariableRef[],GraphConstraintRef[],GraphConstraintRef[])
gethypergraph(graph::AggregateGraph) = graph.hypergraph

mutable struct AggregateNode <: AbstractModelNode
    hypernode::HyperNode
    obj_dict::Dict{Symbol,Any}
    variablemap::Dict{JuMP.VariableRef,JuMP.VariableRef}
    constraintmap::Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}
    nl_constraintmap::Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}
    objective::Union{JuMP.AbstractJuMPScalar,Expr}
end

zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef}))

# StructureGraphs.getstructurenode(node::JuMPNode) = node.structurenode

# #Has constraint references for link constraints
# mutable struct AggragateGraphEdge <: AbstractLinkingEdge
#     structureedge::StructureEdge
#     linkconstraints::Vector{GraphConstraintRef}  #indices in JuMP model of linkconstraints for this edge
# end
StructureGraphs.create_edge(graph::JuMPGraph) = JuMPEdge(StructureEdge(),GraphConstraintRef[])
StructureGraphs.getstructureedge(edge::JuMPEdge) = edge.structureedge

#Construct a structured model, but roll it all into one JuMP model (this is how we solve with JuMP accessible solvers)
function JuMPGraphModel()
    m = JuMP.Model()
    m.ext[:Graph] = JuMPGraph()
    return m
end
is_graphmodel(m::JuMP.Model) = haskey(m.ext,:Graph) ? true : false  #check if the model is a graph model

# Should be defined by base type
# #Add nodes and edges to graph models.  These are used for model instantiation from a graph
function NestedHyperGraphs.add_node!(m::JuMP.Model; index = nv(getgraph(m).hypergraph)+1)
    is_graphmodel(m) || error("Can only add nodes to graph models")
    node = create_node(getgraph(m))
    add_node!(getgraph(m),node,index = index)
    return node
end

function StructureGraphs.add_edge!(m::JuMP.Model,nodes::JuMPNode...)
    is_graphmodel(m) || error("Can only add edges to graph models")
    edge = add_edge!(getgraph(m),nodes...)
    return edge
end

#Define all of the JuMP model extension functions
getgraph(m::JuMP.Model) = haskey(m.ext, :Graph) ? m.ext[:Graph] : error("Model is not a graph model")
StructureGraphs.getnodes(m::JuMP.Model) = getnodes(getgraph(m))
StructureGraphs.getedges(m::JuMP.Model) = getedges(getgraph(m))

StructureGraphs.getnode(m::JuMP.Model,id::Integer) = getnode(getgraph(m),id)  #Grab from the highest level graph if not specified
StructureGraphs.getnode(m::JuMP.Model,sid::Integer,nid::Integer) = getnode(getgraph(m).subgraphlist[sid],nid)

# StructureGraphs.getedge(m::Model,id::LightGraphs.AbstractEdge) = getedge(getgraph(m))[id]
# StructureGraphs.getedge(m::Model,sid::Integer,eid::LightGraphs.AbstractEdge) = getedge(getgraph(m).subgraphlist[sid],eid)

JuMP.objective_function(node::JuMPNode) = node.objective

getnodevariables(node::JuMPNode) = node.variablelist
getnodevariable(node::JuMPNode,index::Integer) = node.variablelist[index]
getnodeconstraints(node::JuMPNode) = node.constraintlist
JuMP.num_variables(node::JuMPNode) = length(node.variablelist)
#get node variables using a symbol lookup
getindex(node::JuMPNode,s::Symbol) = node.obj_dict[s]

#get all of the link constraints from a JuMP model
#Retrieves constraint indices that are link constraints
function getlinkconstraints(m::JuMP.Model)
    is_graphmodel(m) || error("link constraints are only available on graph models")
    return getgraph(m).linkconstraints
end
