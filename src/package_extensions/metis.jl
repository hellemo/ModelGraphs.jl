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

function Metis.partition(ugraph::NodeUnipartiteGraph,proj_map::ProjectionMap,n_parts::Int64; alg = :KWAY, vertex_weights = :model_size, edge_weights = :constraints)

    lightgraph = getlightgraph(ugraph)
    G = Metis.graph(lightgraph)  #Metis Graph object

    #Fill in vertex weights with model size of nodes
    v_weights = C_NULL
    if vertex_weights == :model_size
        v_weights = Vector{Metis.idx_t}(undef, G.nvtxs)
        for (vertex,node) in proj_map.node_map
            v_weights[vertex] = JuMP.num_variables(node)
        end
    end

    #Fill in edge weights with number of linking constraints for each edge
    e_weights = C_NULL
    if edge_weights == :constraints
        e_weights = Vector{Metis.idx_t}(undef, length(G.adjncy))
        adj_start = 1
        edge_count = 0
        for v = 1:G.nvtxs
            adj_end = G.xadj[v+1] - 1
            for out_vertex in G.adjncy[adj_start:adj_end]
                edge = Edge(v,out_vertex)
                edge_count += 1
                if haskey(proj_map.edge_map,edge)
                    hyper_edge = proj_map[edge]
                else
                    edge = Edge(out_vertex,v)
                    hyper_edge = proj_map[edge]
                end
                n_links = sum([length(hyper_edge[i].linkconstraints) for i = 1:length(hyper_edge)])
                e_weights[edge_count] = n_links
            end
        end
    end


    part = Vector{Metis.idx_t}(undef, G.nvtxs)
    edgecut = fill(Metis.idx_t(0), 1)

    if alg === :RECURSIVE
        Metis.METIS_PartGraphRecursive(G.nvtxs, idx_t(1), G.xadj, G.adjncy, v_weights, C_NULL,e_weights,
                                 idx_t(nparts), C_NULL, C_NULL, Metis.options, edgecut, part)
    elseif alg === :KWAY
        Metis.METIS_PartGraphKway(G.nvtxs, Metis.idx_t(1), G.xadj, G.adjncy, v_weights, C_NULL,e_weights,
                            Metis.idx_t(nparts), C_NULL, C_NULL, Metis.options, edgecut, part)
    else
        throw(ArgumentError("unknown algorithm $(repr(alg))"))
    end
    return part

end
