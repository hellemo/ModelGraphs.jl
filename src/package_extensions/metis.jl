import Metis

"""
partition(graph::ModelGraph,n_parts::Int64;alg = :KWAY) --> Vector{Vector{Int64}}

Return a graph partition containing a vector of a vectors of node indices.
"""
# function Metis.partition(graph::ModelGraph,n_parts::Int64,projection::Function = NodeUnipartiteGraph;alg = :KWAY)
#     ugraph = projection(graph)
#
#     lg = getlightgraph(ugraph)
#
#     #TODO Make metis account for weights
#     parts = Metis.partition(lg,n_parts,alg = alg)
#     unique_parts = unique(parts)
#     nparts = length(unique_parts)
#
#     partition_dict = Dict{Int64,Vector{Int64}}((k,[]) for k in unique_parts)
#     for modelnode in getnodes(graph)
#         index = getindex(graph,modelnode)
#         part = parts[index]
#         push!(partition_dict[part],index)
#     end
#
#     partitions = collect(values(partition_dict))
#
#     return partitions
# end

function Metis.partition(ugraph::NodeUnipartiteGraph,nparts::Int64; alg = :KWAY, vertex_weights = :model_size, edge_weights = :n_constraints)

    lightgraph = getlightgraph(ugraph)
    G = Metis.graph(lightgraph)  #Metis Graph object

    #Fill in vertex weights with model size of nodes
    v_weights = C_NULL
    if vertex_weights == :model_size
        v_weights = Vector{Metis.idx_t}(undef, G.nvtxs)
        for (vertex,weight) in ugraph.v_weights
             v_weights[vertex] = weight
        end
    end

    #Fill in edge weights with number of linking constraints for each edge
    e_weights = C_NULL
    if edge_weights == :n_constraints
        e_weights = Vector{Metis.idx_t}(undef, length(G.adjncy))
        adj_start = 1
        edge_count = 0
        for v = 1:G.nvtxs
            v = Int32(v)                #Metis uses 32 bit Integers
            adj_end = G.xadj[v+1] - 1   #This is the positions in G.adjncy corresponding to the current vertex neighbors
            for out_vertex in G.adjncy[adj_start:adj_end]
                edge = Edge(Int64(v),Int64(out_vertex))
                edge_count += 1
                if haskey(ugraph.e_weights,edge)
                    e_weights[edge_count] = ugraph.e_weights[edge]
                else
                    edge = Edge(Int64(out_vertex),Int64(v))
                    @assert haskey(ugraph.e_weights,edge)
                    e_weights[edge_count] = ugraph.e_weights[edge]
                end
            end
            adj_start = adj_end + 1
        end
    end

    part = Vector{Metis.idx_t}(undef, G.nvtxs)
    edgecut = fill(Metis.idx_t(0), 1)

    #Recursive bisection
    if alg === :RECURSIVE
        Metis.METIS_PartGraphRecursive(G.nvtxs, idx_t(1), G.xadj, G.adjncy, v_weights, C_NULL, e_weights, idx_t(nparts), C_NULL, C_NULL, Metis.options, edgecut, part)

    #Multi-level kway partitioning
    elseif alg === :KWAY
        Metis.METIS_PartGraphKway(G.nvtxs, Metis.idx_t(1), G.xadj, G.adjncy, v_weights, C_NULL, e_weights, Metis.idx_t(nparts), C_NULL, C_NULL, Metis.options, edgecut, part)

    else
        throw(ArgumentError("unknown algorithm $(repr(alg))"))
    end
    return part #membership vector
end
