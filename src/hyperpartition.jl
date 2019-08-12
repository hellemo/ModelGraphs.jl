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

struct NodePartition
    nodes::Vector{HyperNode}
    edges::Vector{HyperEdge}
    parent::Union{Nothing,PartitionLayer}
end

struct PartitionParent
    sharednodes::Vector{HyperNode}
    sharededges::Vector{HyperEdge}
    parent::Union{Nothing,PartitionParent}
end

mutable struct HyperPartition
    node_partitions::Vector{NodePartition}  #bottom level partitions
    partition_tree::Vector{PartitionParent}  #tree structure describing recursive structure and shared nodes and edges
end
HyperPartition() = HyperPartition(Vector{NodePartition}(),Vector{PartitionLayer}())

#Convert membership vector to lists of indices
function getpartitionlist(hypergraph::HyperGraph,membership_vector::Vector)
    unique_parts = unique(membership_vector)  #get unique membership entries
    nparts = length(unique_parts)             #number of partitions
    partitions = [Vector{Int64}() for _ = 1:nparts]
    for (vertex,part) in enumerate(membership_vector)
        push!(partitions[part],getnode(hypergraph,vertex))
    end
    return partitions
end

#Naive implementation.  Need to use incidence matrix to do this correctly, but first I need to find better way to deal with subgraph hyperedges
#Make work with subgraphs
function getinducedhyperedges(hypergraph::HyperGraph,nodes::Vector{HyperNode})
    hyperedges = []
    for node in nodes
        for hyperedge in hypergraph.node_map[node]
            if all(hyperedge.vertices in nodes)
                push!(hyperedges,hyperedge)
end

function getcuthyperedges(hypergraph::HyperGraph,partition::Vector)
end


#Simple 2 level partition from a vector of integers
function HyperPartition(hypergraph::AbstractHyperGraph,node_membership_vector::Vector{Int64})
    partition = HyperPartition()

    #convert membership vector to vector of vectors
    hypernode_partitions = getpartitionlist(hypergraph,node_membership_vector)

    #get induced hypergraph from nodes
    partition_edges = getinducedhyperedges.(hypergraph,hypernode_partitions)
    cut_edges = getcuthyperedges(hypergraph,hypernode_partitions)

    node_partitions = Vector{NodePartitions}()
    for i = 1:length(hypernode_partitions)
        push!(node_partitions,NodePartition(hypernode_partitions[i],partition_edges[i]))
    end



    partition.node_partitions = hypernode_partitions
    partitions.


    return hyperpartition

end

#NOTE: Could also be a Dual Clique Graph
function HyperPartition(clique_graph::CliqueExpandedGraph,projection_map::ProjectionMap,membership_vector::Vector{Int64}))

    hyperpartition = HyperPartition()
    #figure out the hypergraph partition based on the graph partition

    return hyperpartition
end

function HyperPartition(bipartite_graph::BipartiteGraph,conversion_map::ProjectionMap,membership_vector::Vector{Int64});selection = :shared_nodes)

    hyperpartition = HyperPartition()

    return hyperpartition
end

function HyperPartition(dual_hyper_graph::AbstractHyperGraph,conversion_map::ProjectionMap,membership_vector::Vector{Int64}))

    hyperpartition = HyperPartition()

    return hyperpartition
end

# #get hypergraphs using induced subgraph
#
# #shared edges cannot be in any partitions
#
# #shared nodes cannot be in any partitions
#
# #shared nodes cannot be incident to a shared edge
#
# #Given a vector of node indices, create a model partition that contains shared edges
# function ModelPartition(hypergraph::HyperGraph,node_membership_vector::Vector{Int64})
#
#     #We need to build a ModelPartition which contains hypergraph partitions, shared edges between partitions, and shared nodes (with their supporting edges)
#     node_partitions,shared_nodes = _identify_partitions(graph)
#
#     return_shared_entities = unique(_map_entities(shared_entities,projection_map))
#     #NOTE: Need to keep vector the same size. #If there are duplicate entries across partitions, then they must also show up in shared
#     return_partition_entities = [unique(_map_entities(local_entitiy,projection_map)) for local_entitiy in local_entities]
#
#     #Make sure no partition entities are in shared entities.  It's possible that a local entity maps to a linkconstraint that is actually shared.
#     for return_part in return_partition_entities
#         filter!(e ->  !(e in return_shared_entities),return_part)
#     end
#
#     #return PartitionData(return_partitions,return_partition_entities,return_shared_entities)
#
#
#     return ModelPartition
# end
