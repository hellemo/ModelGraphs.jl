# #Case 0
# hypergraph = gethypergraph(modelgraph)  #OR getlinkvarhypergraph(modelgraph)  #hypergraph with a node for the master problem.  the linknode gets index 0
# membership_vector = KaHyPar.partition(hypergraph)
# model_partition = ModelPartition(hypergraph,membership_vector)  #create hypergraph partitions, find shared hyperedges
#
# #Case 1
# hypergraph = gethypergraph(model_graph)  #option to add node for master
# clique_graph, conversion_map = clique_expansion(hyper_graph)  #conversion map maps nodes and edges back to hypergraph
# membership_vector = Metis.partition(clique_graph)  #e.g. [1,2,1,1,3,2,2,2,3,4,4,4]
# model_partition = ModelPartition(clique_graph,conversion_map,membership_vector)
#
# #Case 2
# hypergraph = gethypergraph(model_graph)
# dual_graph, conversion_map = dual_hyper_graph(hypergraph)
# membership_vector = KaHyPar.partition(dual_graph,4)
# model_partition = ModelPartition(dual_hyper_graph,conversion_map,membership_vector)
#
# #Case 3
# hypergraph = gethypergraph(model_graph)
# bipartite_graph, conversion_map = star_expansion(hypergraph)
# membership_vector = Metis.partition(bipartite_graph,4)
# model_partition = ModelPartition(bipartite_graph,conversion_map,membership_vector;selection = :shared_nodes)
#
# #Case 4
# hypergraph = gethypergraph(model_graph)
# dual_clique_graph, conversion_map = dual_clique_expansion(hypergraph)
# membership_vector = Metis.partition(dual_clique_graph,4)
# model_partition = ModelPartition(dual_clique_graph,conversion_map,membership_vector)
mutable struct ModelPartition
    partitions::Vector{HyperGraph}
    sharednodes::Union{Nothing,Vector{HyperNode}}
    sharededges::Union{Nothing,Vector{HyperEdge}}
    subpartitions::Union{Nothing,Vector{ModelPartition}}  #1 -> part1, 2-> part2
end

#Constructors
function ModelPartition(hypergraph::HyperGraph,node_membership_vector::Vector{Int64})
    modelpartition = ModelPartition()
    #get hypernodes

    #get partition of nodes

    #we need to create the vertex induced subgraph given each node partition

    #shared edges cannot be in any partitions

    #shared nodes cannot be in any partitions

    #do the same thing for subpartitions


end

#Given a vector of node indices, create a model partition that contains shared edges
function ModelPartition(hypergraph::HyperGraph,node_membership_vector::Vector{Int64})

    #We need to build a ModelPartition which contains hypergraph partitions, shared edges between partitions, and shared nodes (with their supporting edges)
    node_partitions,shared_nodes = _identify_partitions(graph)

    return_shared_entities = unique(_map_entities(shared_entities,projection_map))
    #NOTE: Need to keep vector the same size. #If there are duplicate entries across partitions, then they must also show up in shared
    return_partition_entities = [unique(_map_entities(local_entitiy,projection_map)) for local_entitiy in local_entities]

    #Make sure no partition entities are in shared entities.  It's possible that a local entity maps to a linkconstraint that is actually shared.
    for return_part in return_partition_entities
        filter!(e ->  !(e in return_shared_entities),return_part)
    end

    return PartitionData(return_partitions,return_partition_entities,return_shared_entities)


    return ModelPartition
end


#Given a vector of node indices, create a model partition that contains shared edges
#Hyperedges are numbered in a graph
function ModelPartition(graph::HyperGraph,conversion_map::ConversionMap)

    #Input checking

    shared_nodes = nothing

    #create hypergraph object for each partition


    #find shared edges


    #create subpartitions

end
