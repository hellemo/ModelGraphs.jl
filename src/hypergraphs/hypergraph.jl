import Base: ==,string,print,show
using LightGraphs

abstract type AbstractHyperGraph <: LightGraphs.AbstractGraph{Int64} end
abstract type AbstractHyperEdge <: LightGraphs.AbstractEdge{Int64} end

struct HyperNode
    index::Int64
end
# HyperNode(index::Int64) = HyperNode()

struct HyperEdge <: AbstractHyperEdge
    index::Int64
    vertices::Set{HyperNode}
end
#HyperEdge(index::Int64,t::Set{HyperNode}) = HyperEdge(Dict{AbstractHyperGraph,Int64}(hypergraph => index),t)
HyperEdge(index::Int64,t::Vector{HyperNode}) = HyperEdge(index,Set(t))
HyperEdge(index::Int64,t::HyperNode...) = HyperEdge(index,Set(collect(t)))

# LightGraphs.AbstractGraph
mutable struct HyperGraph <: AbstractHyperGraph
    vertices::Vector{HyperNode}
    hyperedge_map::OrderedDict{Int64,HyperEdge}                  #look up hyperedges by index in the hypergraph
    hyperedges::OrderedDict{Set,HyperEdge}       #look up hyperedges using hypernodes.  These are LOCAL to the hypergraph
    node_map::Dict{HyperNode,Vector{HyperEdge}}  #map hypernodes to hyperedges they are incident to
end
HyperGraph() = HyperGraph(HyperNode[],OrderedDict{Int64,HyperEdge}(),OrderedDict{Set,AbstractHyperEdge}(),Dict{HyperNode,Vector{AbstractHyperEdge}}())

#SparseMatrix from hypergraph
function SparseArrays.sparse(hypergraph::HyperGraph)
    #Build up I and J.  Assume V = 1
    I = []
    J = []
    for hyperedge in gethyperedges(hypergraph)
        edge_index = getindex(hypergraph,hyperedge)
        node_indices = []
        #NOTE: Order shouldn't matter here, but need to check
        for hypernode in hyperedge.vertices
            node_index = getindex(hypergraph,hypernode)
            push!(node_indices,node_index)
        end
        for node_index in sort(node_indices)
            push!(I,node_index)
            push!(J,edge_index)
        end
    end
    V = Int.(ones(length(I)))
    return SparseArrays.sparse(I,J,V)
end

#HyperNode
function LightGraphs.add_vertex!(hypergraph::AbstractHyperGraph)
    (nv(hypergraph) + one(Int) <= nv(hypergraph)) && return false       # test for overflow
    v = length(hypergraph.vertices)+1
    hypernode = HyperNode(v)
    push!(hypergraph.vertices,hypernode)
    hypergraph.node_map[hypernode] = HyperEdge[]
    return hypernode
end

add_node!(hypergraph::AbstractHyperGraph) = LightGraphs.add_vertex!(hypergraph)


getnode(hypergraph::HyperGraph,index::Int64) = hypergraph.vertices[index]
Base.getindex(hypergraph::HyperGraph,node::HyperNode) = node.index#_map[hypergraph]
getnodes(hypergraph::HyperGraph) = hypergraph.vertices


#HyperEdge
Base.reverse(e::HyperEdge) = "hyperedge doesn't not support reverse"
==(h1::HyperEdge,h2::HyperEdge) = collect(h1.vertices) ==  collect(h2.vertices)  #vertices are sorted when added

LightGraphs.add_edge!(graph::AbstractHyperGraph,vertices::Int...) = add_hyperedge!(graph,vertices...)

gethypernodes(edge::HyperEdge) = collect(edge.vertices)

#Add new LOCAL HyperEdge to a HyperGraph
function add_hyperedge!(hypergraph::AbstractHyperGraph,vertices::Int64...)
    hypernodes = map(x -> getnode(hypergraph,x),vertices)
    return add_hyperedge!(hypergraph,hypernodes...)
end

function add_hyperedge!(hypergraph::AbstractHyperGraph,hypernodes::HyperNode...)
    @assert length(hypernodes) > 1
    hypernodes = Set(collect(hypernodes))
    if has_edge(hypergraph,hypernodes)
        return hypergraph.hyperedges[hypernodes]
    else
        index = ne(hypergraph) + 1
        hyperedge = HyperEdge(index,hypernodes...)
        for hypernode in hypernodes
            push!(hypergraph.node_map[hypernode], hyperedge)
        end
        hypergraph.hyperedges[hypernodes] = hyperedge
        hypergraph.hyperedge_map[index] = hyperedge
        return hyperedge
    end
end
#Get hyperedges

function gethyperedge(hypergraph::AbstractHyperGraph,vertices::Int64...)
    hypernodes = map(x -> getnode(hypergraph,x),vertices)
    return hypergraph.hyperedges[Set(hypernodes)]
end
gethyperedge(hypergraph::HyperGraph,edge_index::Int64) = hypergraph.hyperedge_map[edge_index]
gethyperedges(hypergraph::AbstractHyperGraph) = values(hypergraph.hyperedges)   #local hyperedges
#getallhyperedges(hypergraph::AbstractHyperGraph) = values(hypergraph.hyperedge_map)
getedges(hypergraph::AbstractHyperGraph) = gethyperedges(hypergraph)
vertices(hyperedge::HyperEdge) = collect(hyperedge.vertices)
Base.getindex(hypergraph::HyperGraph,edge::HyperEdge) = edge.index#_map[hypergraph]

#LightGraphs Interface
LightGraphs.edges(graph::AbstractHyperGraph) = graph.hyperedges #HyperEdgeIter(g)
LightGraphs.edgetype(graph::AbstractHyperGraph) = HyperEdge
LightGraphs.has_edge(graph::AbstractHyperGraph,edge::HyperEdge) = edge in values(graph.hyperedges)
function LightGraphs.has_edge(graph::AbstractHyperGraph,hypernodes::Set{HyperNode})
    return haskey(graph.hyperedges,hypernodes)
end

LightGraphs.has_vertex(graph::AbstractHyperGraph, v::Integer) = v in vertices(graph)
LightGraphs.is_directed(graph::AbstractHyperGraph) = false
LightGraphs.is_directed(::Type{AbstractHyperGraph}) = false
LightGraphs.ne(graph::AbstractHyperGraph) = length(graph.hyperedge_map)
LightGraphs.nv(graph::AbstractHyperGraph) = length(graph.vertices)
LightGraphs.vertices(graph::AbstractHyperGraph) = graph.vertices

# LightGraphs.rem_edge!
#TODO This shouldn't be too bad
LightGraphs.rem_edge!(g::AbstractHyperGraph,e::HyperEdge) = throw(error("Edge removal not yet supported on hypergraphs"))

# LightGraphs.rem_vertex!
#TODO Delete any associated edges with the vertex
LightGraphs.rem_vertex!(g::AbstractHyperGraph) = throw(error("Vertex removal not yet supported on hypergraphs"))

#NOTE Inefficient neighbors implementation
#Could use a SparseArray to do this faster
function LightGraphs.all_neighbors(g::HyperGraph,node::HyperNode)
    hyperedges = g.node_map[node]  #incident hyperedges to the hypernode
    neighbors = []
    for edge in hyperedges
        append!(neighbors,[vert for vert in edge.vertices if vert != node])
    end
    return unique(neighbors)
end

function incident_edges(g::HyperGraph,node::HyperNode)
    hyperedges = g.node_map[node]  #incident hyperedges to the hypernode
    return hyperedges
end

#TODO Get neighborhood of a node out to a given distance
function neighborhood(g::HyperGraph,node::HyperNode,distance::Int64)
    dist = 0
    neighbor_list = []
    edge_list = []

    nodes_to_check = [node]
    checked_nodes = []
    checked_edges = []
    while dist < distance
        if !isempty(nodes_to_check)
            for node in nodes_to_check
                neighbors = LightGraphs.all_neighbors(g,node)
                iedges = incident_edges[g,node]

                f1 = (node) -> !(node in checked_nodes)
                f2 = (edge) -> !(edge in checked_edges)

                filter!(f1,neighbors)  #these are the new neighbors
                filter!(f2,iedges)

                popfirst!(nodes_to_check)
                append!(nodes_to_check,neighbors)
                push!(checked_nodes,node)

            end
            dist += 1
        end
    end

    return nodes,edges
end


####################################
#Print Functions
####################################
function string(graph::AbstractHyperGraph)
    "Hypergraph: "*"($(nv(graph)) , $(ne(graph)))"
end
print(io::IO, graph::AbstractHyperGraph) = print(io, string(graph))
show(io::IO,graph::AbstractHyperGraph) = print(io,graph)

function string(node::HyperNode)
    "Hypernode: "*"$(node.index)" #_map))"
end
print(io::IO,node::HyperNode) = print(io, string(node))
show(io::IO,node::HyperNode) = print(io,node)

function string(edge::HyperEdge)
    "HyperEdge: "*"$(collect(edge.vertices))"
end
print(io::IO,edge::HyperEdge) = print(io, string(edge))
show(io::IO,edge::HyperEdge) = print(io,edge)



# #Copy hypergraph.  Retain subgraphs too
# function Base.copy(hypergraph::HyperGraph)
#     copy_hypergraph = HyperGraph()
#     for node in getnodes(hypergraph)
#         add_node!(copy_hypergraph)
#     end
#     for edge in get
# end





#LightGraphs.degree(g::HyperGraph,v::Int) = length(all_neighbors(g,v))

# function rem_edge!(g::SimpleGraph, e::SimpleGraphEdge)
#     i = searchsorted(g.fadjlist[src(e)], dst(e))
#     isempty(i) && return false   # edge not in graph
#     j = first(i)
#     deleteat!(g.fadjlist[src(e)], j)
#     if src(e) != dst(e)     # not a self loop
#         j = searchsortedfirst(g.fadjlist[dst(e)], src(e))
#         deleteat!(g.fadjlist[dst(e)], j)
#     end
#     g.ne -= 1
#     return true # edge successfully removed
# end

# function LightGraphs.add_edge!(g::HyperGraph,edge::SimpleEdge)
#     hedge = HyperEdge(src(edge),dst(edge))
#     inserted = add_edge!(g,hedge)
#     return inserted
# end

# h_index = 0
# for hyperedge in gethyperedges(hypergraph)
#     h_index += 1
#     for hypernode in hyperedge.vertices
#         vertex = getindex(hypergraph,hypernode)
#         push!(I,vertex)
#         push!(J,h_index)
#     end
# end
#
# #h_index = ne(hypergraph) #need to track hyperedges in higher level graph
# for subgraph in subgraphs(hypergraph)
#     for hyperedge in gethyperedges(subgraph)
#         h_index += 1
#         for hypernode in hyperedge.vertices
#             hypergraph_vertex = getindex(hypergraph,hypernode)
#             push!(I,hypergraph_vertex)
#             push!(J,h_index)
#         end
#     end
# end

# function add_node!(hypergraph::AbstractHyperGraph,node::HyperNode)
#     push!(hypergraph.vertices,node)
#     v = length(hypergraph.vertices)
#     node.index = v
#     hypergraph.node_map[node] = HyperEdge[]
#     return node
# end

# #Subgraph
# subgraphs(graph::HyperGraph) = graph.subgraphs
# getsubgraphs(graph::HyperGraph) = subgraphs(graph)
#
# function add_subgraph!(graph::HyperGraph,subgraph::HyperGraph)
#     push!(graph.subgraphs,subgraph)
#     subgraph.index = length(subgraphs(graph))           #subgraph index
#     for node in getnodes(subgraph)
#         add_node!(graph,node)  #add the HyperNode to the hypergraph
#     end
#
#     #TODO Fix.  This is slow.
#     for edge in getedges(subgraph)
#         add_sub_hyperedge!(graph,edge)
#     end
#     for subgraph in subgraphs(subgraph)
#         add_subgraph!(graph,subgraph)
#     end
#     return graph
# end
