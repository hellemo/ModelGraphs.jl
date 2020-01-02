
#Create a hypergraph object from a ModelGraph.  Return the hypergraph and a mapping of hypergraph nodes and edges to modelgraph modelnodes and linkedges.
#A hypergraph has topology functions and partitioning interfaces.
function gethypergraph(graph::ModelGraph;include_master_nodes = false)
    hypergraph = HyperGraph()

    hyper_map = Dict()  #mapping from hypergraph nodes to modelnodes and link_edges
    node_map = Dict()

    for node in all_nodes(graph)
        hypernode = add_node!(hypergraph)
        hyper_map[hypernode] = node
        node_map[node] = hypernode
    end

    for edge in all_edges(graph)
        nodes = edge.nodes
        hypernodes = [node_map[modelnode] for modelnode in nodes]
        hyperedge = add_hyperedge!(hypergraph,hypernodes...)
        hyper_map[hyperedge] = edge
    end

    #Add master nodes to the hypergraph
    if include_master_nodes == true
        nothing
    end


    return hypergraph,hyper_map
end
