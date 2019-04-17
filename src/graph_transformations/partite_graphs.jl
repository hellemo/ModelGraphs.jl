#A Bipartite graph (Standard Graph) where one set of nodes corresponds to
mutable struct BipartiteGraph <: AbstractModelGraph
    structuregraph::StructureGraphs.StructureGraph
    part1::Vector{Int64}   #partition 1 of the bipartite graph
    part2::Vector{Int64}   #partition 2 of the bipartite graph
end

BipartiteGraph() = BipartiteGraph(StructureGraph(LightGraphs.Graph()),Int64[],Int64[])

StructureGraphs.create_node(graph::BipartiteGraph) = StructureNode()
StructureGraphs.create_edge(graph::BipartiteGraph) = StructureEdge()

function string(graph::BipartiteGraph)
    "Bipartite Graph\ngraph_id: "*string(getlabel(graph))*"\nnodes:"*string((length(getnodes(graph))))
end

#A Unipartite graph (Standard Graph) where nodes correspond to model nodes and edges correspond to links between model nodes.
#This structure can be used to convert a ModelGraph structure to a Shared Constraint structure
mutable struct NodeUnipartiteGraph <: AbstractModelGraph
    structuregraph::StructureGraphs.StructureGraph
    v_weights::Dict{Int64,Int64}                     #vertex weights
    e_weights::Dict{LightGraphs.AbstractEdge,Int64}  #edge weights
end

NodeUnipartiteGraph() = NodeUnipartiteGraph(StructureGraph(LightGraphs.Graph),Dict{StructureNode,Int64}(),Dict{StructureNode,Int64}())

StructureGraphs.create_node(graph::NodeUnipartiteGraph) = StructureNode()
StructureGraphs.create_edge(graph::NodeUnipartiteGraph) = StructureEdge()

function string(graph::NodeUnipartiteGraph)
    "Unipartite Graph\ngraph_id: "*string(getlabel(graph))*"\nnodes:"*string((length(getnodes(graph))))
end

#A Unipartite graph (Standard Graph) where nodes correspond to link constraints and edges correspond to shared model nodes.
#This structure can be used to convert a ModelGraph structure to a Shared Model (Variable) structure
mutable struct LinkUnipartiteGraph <: AbstractModelGraph
    structuregraph::StructureGraphs.StructureGraph
    v_weights::Dict{Int64,Int64}                     #vertex weights
    e_weights::Dict{LightGraphs.AbstractEdge,Int64}  #edge weights
end

LinkUnipartiteGraph() = LinkUnipartiteGraph(StructureGraph(LightGraphs.Graph),Dict{StructureNode,Int64}(),Dict{StructureNode,Int64}())

StructureGraphs.create_node(graph::LinkUnipartiteGraph) = StructureNode()
StructureGraphs.create_edge(graph::LinkUnipartiteGraph) = StructureEdge()

function string(graph::LinkUnipartiteGraph)
    "Unipartite Graph\ngraph_id: "*string(getlabel(graph))*"\nnodes:"*string((length(getnodes(graph))))
end
