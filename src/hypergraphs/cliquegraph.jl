import LightGraphs.SimpleGraphs.SimpleEdge
#A simple graph that keeps track of how many duplicate edges have been added
mutable struct CliqueExpandedGraph <: LightGraphs.AbstractGraph{Int64}
    lightgraph::LightGraphs.Graph
    edge_counts::Dict{SimpleEdge,Int64}
end
CliqueExpandedGraph() = CliqueExpandedGraph(LightGraphs.Graph(),Dict{SimpleEdge,Int64}())


#TODO: Replace with simple macro
LightGraphs.add_vertex!(graph::CliqueExpandedGraph) = LightGraphs.add_vertex!(graph.lightgraph)


function LightGraphs.add_edge!(graph::CliqueExpandedGraph,from::Int64,to::Int64)
    if has_edge(graph,from,to)
        graph.edge_counts[SimpleEdge(sort([from,to])...)] += 1
        inserted = false
    else
        LightGraphs.add_edge!(graph.lightgraph,from,to)
        graph.edge_counts[SimpleEdge(sort([from,to])...)] = 1
        inserted = true
    end
    return inserted
end


LightGraphs.edges(graph::CliqueExpandedGraph) = LightGraph.edges(graph.lightgraph)
LightGraphs.edgetype(graph::CliqueExpandedGraph) = LightGraphs.SimpleGraphs.SimpleEdge{Int64}
LightGraphs.has_edge(graph::CliqueExpandedGraph,from::Int64,to::Int64) = LightGraphs.has_edge(graph.lightgraph,from,to)


LightGraphs.has_vertex(graph::CliqueExpandedGraph, v::Integer) = LightGraphs.has_vertex(graph.lightgraph,v)
LightGraphs.is_directed(graph::CliqueExpandedGraph) = false
LightGraphs.is_directed(::Type{CliqueExpandedGraph}) = false
LightGraphs.ne(graph::CliqueExpandedGraph) = LightGraphs.ne(graph.lightgraph)
LightGraphs.nv(graph::CliqueExpandedGraph) = LightGraphs.nv(graph.lightgraph)
LightGraphs.vertices(graph::CliqueExpandedGraph) = LightGraphs.vertices(graph.lightgraph)
