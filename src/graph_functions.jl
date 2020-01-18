#Retrieve the neighbors outside of a subgraph
function neighborhood_expansion(mg::ModelGraph,subgraph::ModelGraph,boundary_nodes::Vector{ModelNodes},overlap::Int64)
    hypergraph,hypermap = gethypergraph(mg)

    modelnodes= = ModelNode[]
    subgraph_nodes = [hypermap[node] for node in subgraph.modelnodes]  #these are hypernodes

    for node in boundary_nodes
        hypernode = hypermap[node]
        neighbors = neighborhood(hypergraph,hypernode)
    end


    return modelnodes,linkedges,boundary_edges
end
