function LightGraphs.all_neighbors(mg::ModelGraph,node::ModelNode)
    hypergraph,hypermap = gethypergraph(mg)
    hypernode = hypermap[node]
    neighbors = LightGraphs.all_neighbors(hypergraph,hypernode)
    return [hypermap[neighbor] for neighbor in neighbors]
end

function incident_edges(mg::ModelGraph,node::ModelNode)
    hypergraph,hypermap = gethypergraph(mg)
    hypernode = hypermap[node]
    i_edges = incident_edges(hypergraph,hypernode)
    return [hypermap[edge] for edge in i_edges]
end

function neighborhood(mg::ModelGraph,node::ModelNode,distance::Int64;exclude_nodes = ModelNode[])
    #Map to hypergraph
    hypergraph,hypermap = gethypergraph(mg)
    hypernode = hypermap[node]
    exclude_hypernodes = [hypermap[node] for node in exclude_nodes]

    #Get neighbors.  We get a list for each layer of distance
    neighbor_list,edge_list = neighborhood(hypergraph,hypernode,distance;exclude_nodes = exclude_hypernodes)

    #Map back to ModelGraph
    f = (node) -> hypermap[node]
    return_neighbor_list = [f.(l) for l in neighbor_list]
    return_edge_list = [f.(l) for l in edge_list]
    return return_neighbor_list,return_edge_list
end

#Perform neighborhood expansion on a subgraph boundary
function neighborhood_expansion(mg::ModelGraph,subgraph::ModelGraph,boundary_nodes::Vector{ModelNode},overlap::Int64)
    #hypergraph,hypermap = gethypergraph(mg)
    #subgraph_nodes = [hypermap[node] for node in subgraph.modelnodes]  #these are subgraph hypernodes
    subgraph_nodes = subgraph.modelnodes
    new_modelnodes = ModelNode[]
    new_linkedges = LinkEdge[]
    boundary_edges = LinkEdge[]

    for node in boundary_nodes
        # neighbor_node_list,edge_list = neighborhood(hypergraph,hypernode,overlap;exclude = subgraph_nodes)
        neighbor_node_list,edge_list = neighborhood(mg,node,overlap,exclude_nodes = subgraph_nodes)

        append!(new_modelnodes,vcat(neighbor_node_list...))
        append!(new_linkedges,vcat(edge_list...))
        #now get the incident edges from the last layer of neighbor_nodes.  These are the link constraints we need to coordinate.

        final_nodes = neighbor_node_list[end]
        new_i_edges = []
        for node in final_nodes
            i_edges = incident_edges(mg,node)
            f = (edge) -> !(edge in new_linkedges)
            filter!(f,i_edges)  #don't include edges that are internal to the subproblem.  We just want the boundary.
            append!(new_i_edges,i_edges)
        end
        append!(boundary_edges,new_i_edges)
    end

    #Now just need final boundary edges from outer layers
    return new_modelnodes, new_linkedges, boundary_edges
end
