function gethypergraph(graph::ModelGraph)
    hypergraph = HyperGraph()
    hyper_map = Dict()  #mapping from hypergraph nodes to modelnodes and link_edges
    for node in all_nodes(graph)
        hypernode = add_vertex!(hypergraph)
        hyper_map[hypernode] = node
    end

    for (i,edge) in enumerate(all_edges(graph))
        nodes = edge.nodes
        hypernodes = [hyper_map[node] for node in nodes]
    end

    return hypergraph,hyper_map
end
