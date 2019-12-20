function gethypergraph(graph::ModelGraph)
    hypergraph = HyperGraph()
    hyper_map = Dict()  #mapping from hypergraph nodes to modelnodes and link_edges
    for node in getnodes(graph)
        hypernode = add_vertex!(hypergraph)
    end



    return hypergraph,hyper_map
end
