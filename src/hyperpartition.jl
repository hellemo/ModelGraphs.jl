abstract type AbstractPartition end

struct SubgraphPartition <: AbstractPartition
    hypergraph::HyperGraph
    parent::Union{Nothing,AbstractPartition}
end

struct PartitionParent <: AbstractPartition
    sharednodes::Vector{HyperNode}          #think master node with link variables
    sharededges::Vector{HyperEdge}          #link constraints
    parent::Union{Nothing,PartitionParent}
    children::Vector{AbstractPartition}
end
PartitionParent(sharednodes::Vector{HyperNode},sharededges::Vector{HyperEdge}) = PartitionParent(sharednodes,sharededges,nothing,Vector{AbstractPartition}())
PartitionParent(sharededges::Vector{HyperEdge}) = PartitionParent(Vector{HyperNode}(),sharededges,nothing,Vector{AbstractPartition}())

mutable struct HyperPartition
    partitions::Vector{SubgraphPartition}  #bottom level partitions
    parents::Vector{PartitionParent}  #tree structure describing recursive structure and shared nodes and edges
end
HyperPartition() = HyperPartition(Vector{SubgraphPartition}(),Vector{PartitionParent}())

#Convert membership vector to lists of indices
function getpartitionlist(hypergraph::HyperGraph,membership_vector::Vector)
    unique_parts = unique(membership_vector)  #get unique membership entries

    nparts = length(unique_parts)             #number of partitions

    #partitions = [Vector{HyperNode}() for _ = 1:nparts]
    partitions = Dict{Int64,Vector{HyperNode}}((k,[]) for k in unique_parts)
    for (vertex,part) in enumerate(membership_vector)
        push!(partitions[part],getnode(hypergraph,vertex))
    end

    return collect(values(partitions))
end

#NOTE Naive implementation to get induced and shared hyperedges given a set of node partitions
#TODO: Use incidence matrix and partition list to figure out what the induced and cut edges are
function identifyhyperedges(hypergraph::HyperGraph,partitions::Vector{Vector{HyperNode}})
    nparts = length(partitions)
    induced_edges = [Vector{HyperEdge}() for _ = 1:nparts]
    shared_edges = Vector{HyperEdge}()  #between partitions

    checked_edges = Vector{HyperEdge}()

    for (i,partition) in enumerate(partitions)
        for hypernode in partition
            for hyperedge in hypergraph.node_map[hypernode]  #getedges(hypergraph,hypernode)
                if !(hyperedge in checked_edges)  #If it's a new link constraint
                    edge_hypernodes = collect(hyperedge.vertices)
                    if all(node -> node in partition,edge_hypernodes)
                        push!(induced_edges[i],hyperedge)
                    else
                        push!(shared_edges,hyperedge)
                    end
                    push!(checked_edges,hyperedge)
                end
            end

            for subgraph in subgraphs(hypergraph)
                if haskey(subgraph.node_map,hypernode)
                    for hyperedge in subgraph.node_map[hypernode]
                        if !(hyperedge in checked_edges)  #If it's a new link constraint
                            edge_hypernodes = collect(hyperedge.vertices)
                            if all(node -> node in partition,edge_hypernodes)
                                push!(induced_edges[i],hyperedge)
                            else
                                push!(shared_edges,hyperedge)
                            end
                            push!(checked_edges,hyperedge)
                        end
                    end
                end
            end
        end
    end
    return induced_edges,shared_edges
end

#Simple 2 level partition from a vector of integers
function HyperPartition(hypergraph::NHG.AbstractHyperGraph,node_membership_vector::Vector{Int64})
    hyperpartition = HyperPartition()

    #convert membership vector to vector of vectors
    hypernode_vectors = getpartitionlist(hypergraph,node_membership_vector)
    induced_edge_partitions,shared_edges = identifyhyperedges(hypergraph,hypernode_vectors)

    #Create new Hypergraphs
    new_hypers = Vector{HyperGraph}()
    for i = 1:length(hypernode_vectors)
        hypernodes = hypernode_vectors[i]
        induced_edges = induced_edge_partitions[i]
        hyper = HyperGraph()

        for hypernode in hypernodes
            add_node!(hyper,hypernode)
        end

        for hyperedge in induced_edges
            add_hyperedge!(hyper,hyperedge)
        end

        push!(new_hypers,hyper)
    end

    partition_parent = PartitionParent(shared_edges)
    partitions = Vector{SubgraphPartition}()
    for i = 1:length(new_hypers)
        push!(partitions,SubgraphPartition(new_hypers[i],partition_parent))
    end
    hyperpartition.partitions = partitions
    hyperpartition.parents = [partition_parent]

    return hyperpartition
end

#Create a SubModelGraph from a HyperGraph specification
function create_sub_modelgraph(modelgraph::ModelGraph,hypergraph::HyperGraph)
    submg = ModelGraph()
    submg.hypergraph = hypergraph

    for hypernode in getnodes(hypergraph)
        modelnode = getnode(modelgraph,hypernode)
        submg.modelnodes[hypernode] = modelnode
    end

    i = 1
    for hyperedge in getedges(hypergraph)
        linkedge = findlinkedge(modelgraph,hyperedge)  #could be in a subgraph
        submg.linkedges[hyperedge] = linkedge
        for linkconstraintref in linkedge.linkconstraints
            linkconstraint = LinkConstraint(linkconstraintref)
            #idx = linkconstraintref.idx
            submg.linkconstraints[i] = linkconstraint
            i += 1
        end
    end
    return submg
end


#TODO
# function HyperPartition(clique_graph::NHG.CliqueExpandedGraph,projection_map::NHG.ProjectionMap,membership_vector::Vector{Int64}) #NOTE: Could also be a Dual Clique Graph
#
#     hyperpartition = HyperPartition()
#
#     #figure out the hypergraph partition based on the graph partition
#
#     return hyperpartition
# end
#
# function HyperPartition(bipartite_graph::NHG.BipartiteGraph,projection_map::NHG.ProjectionMap,membership_vector::Vector{Int64};selection = :shared_nodes)
#
#     hyperpartition = HyperPartition()
#
#     return hyperpartition
# end
#
# function HyperPartition(dual_hyper_graph::AbstractHyperGraph,projection_map::NHG.ProjectionMap,membership_vector::Vector{Int64})
#
#     hyperpartition = HyperPartition()
#
#     return hyperpartition
# end
#















































# #Naive implementation.  Need to use incidence matrix to do this correctly, but first I need to find better way to deal with subgraph hyperedges
#get induced hypergraph from nodes
# induced_edges = getinducedhyperedges.(hypergraph,hypernode_partitions)
# cut_edges = getcuthyperedges(hypergraph,hypernode_partitions)
# function getinducedhyperedges(hypergraph::HyperGraph,nodes::Vector{HyperNode})
#     hyperedges = []
#     for node in nodes
#         for hyperedge in hypergraph.node_map[node]  #this is only edges in the given hypergraph
#             if all(hyperedge.vertices in nodes)
#                 push!(hyperedges,hyperedge)
#             end
#         end
#     end
#     return hyperedges
# end
#
# function getcuthyperedges(hypergraph::HyperGraph,partition::Vector{Vector{HyperNode}})
#
# end
# #Case 0
# hypergraph = gethypergraph(modelgraph)  #OR getlinkvarhypergraph(modelgraph)  #hypergraph with a node for the master problem.  the linknode gets index 0
# membership_vector = KaHyPar.partition(hypergraph)
# model_partition = ModelPartition(hypergraph,membership_vector)  #create hypergraph partitions, find shared hyperedges
#
# #Case 1
# hypergraph = gethypergraph(model_graph)  #option to add node for master
# clique_graph, projection_map = clique_expansion(hyper_graph)  #conversion map maps nodes and edges back to hypergraph
# membership_vector = Metis.partition(clique_graph)  #e.g. [1,2,1,1,3,2,2,2,3,4,4,4]
# hyper_partition = HyperPartition(clique_graph,projection_map,membership_vector)
# new_graph, agg_map = aggregate(modelgraph,hyper_partition)
#
# #Case 2
# hypergraph = gethypergraph(model_graph)
# dual_graph, projection_map = dual_hyper_graph(hypergraph)
# membership_vector = KaHyPar.partition(dual_graph,4)
# model_partition = ModelPartition(dual_hyper_graph,projection_map,membership_vector)
#
# #Case 3
# hypergraph = gethypergraph(model_graph)
# bipartite_graph, projection_map = star_expansion(hypergraph)
# membership_vector = Metis.partition(bipartite_graph,4)
# model_partition = ModelPartition(bipartite_graph,projection_map,membership_vector;selection = :shared_nodes)
#
# #Case 4
# hypergraph = gethypergraph(model_graph)
# dual_clique_graph, projection_map = dual_clique_expansion(hypergraph)
# membership_vector = Metis.partition(dual_clique_graph,4)
# model_partition = ModelPartition(dual_clique_graph,projection_map,membership_vector)
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
