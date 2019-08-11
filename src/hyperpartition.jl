# #Case 0
# hypergraph = gethypergraph(modelgraph)  #OR getlinkvarhypergraph(modelgraph)  #hypergraph with a node for the master problem.  the linknode gets index 0
# membership_vector = KaHyPar.partition(hypergraph)
# model_partition = ModelPartition(hypergraph,membership_vector)  #create hypergraph partitions, find shared hyperedges
#
# #Case 1
# hypergraph = gethypergraph(model_graph)  #option to add node for master
# clique_graph, conversion_map = clique_expansion(hyper_graph)  #conversion map maps nodes and edges back to hypergraph
# membership_vector = Metis.partition(clique_graph)  #e.g. [1,2,1,1,3,2,2,2,3,4,4,4]
# hyper_partition = HyperPartition(clique_graph,conversion_map,membership_vector)
# new_graph, agg_map = aggregate(modelgraph,hyper_partition)
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

# mutable struct ModelPartition
#     partitions::Vector{Vector}          #use induced subgraphs to get the internal edges
#     sharednodes::Vector                 #shared between partitions (e.g. vertex separators)
#     sharededges::Vector                 #span partitions (not connected to a shared node)
# end

#subpartitions::Union{Nothing,Vector{ModelPartition}}  #1 -> part1, 2-> part2


# partitions = [[1,2],[3,4],[5,6],[7,8]]  #subpartitions
# sharednodes = [9, [10],[11],[12],[13]]  #shared nodes at each level
# sharededges = [1, [2],[3],[4],[5]]
#
# partitions = [[[1,2,3],[4,5,6]],[7,8,9],[10,11,12],[13,14,15]]]
# sharedhyperodes = [9, [10],[11],[12],[13]]
# sharedhyperedges = [1, [2],[3],[4],[5]]

struct NodePartition
    nodes::Vector{AbstractHyperGraph}
    parent::PartitionLayer
end

#TODO
function getpartitions(hypergraph::AbstractHyperGraph,node_membership_vector::Vector{Int64})
end

struct SharedEntities
    sharednodes::Vector{HyperNode}
    sharededges::Vector{HyperEdge}
    parent::Union{Nothing,PartitionLayer}
end

mutable struct HyperPartition
    node_partitions::Vector{NodePartition}  #bottom level partitions
    partition_tree::Vector{PartitionLayer}  #tree structure describing recursive structure and shared nodes and edges
end

#Simple 2 level partition from a vector of integers
function HyperPartition(hypergraph::AbstractHyperGraph,node_membership_vector::Vector{Int64})
    hyperpartition = HyperPartition()

    #convert membership vector to vector of vectors
    hypernode_partitions = getpartitions(hypergraph,node_membership_vector)  #Vector of NodePartition

    
    hyperpartition.hypergraph_paritions = getinducedsubgraph.(hypernode_partitions)

    #We have to do this if edge cuts are not given
    hyperpartition.shared_edges = find_edge_cuts(hypergraph,hypernode_partitions)

    return hyperpartition

end

function HyperPartition(clique_graph::CliqueExpandedGraph,conversion_map::ProjectionMap,membership_vector::Vector{Int64}))

    hyperpartition = HyperPartition()

    #figure out the hypergraph partition based on the graph partition




    return hyperpartition
end

#get hypergraphs using induced subgraph

#shared edges cannot be in any partitions

#shared nodes cannot be in any partitions

#shared nodes cannot be incident to a shared edge

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

    #return PartitionData(return_partitions,return_partition_entities,return_shared_entities)


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
