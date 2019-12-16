#NOTE The Partition object describes recursive hypergraph partitions.  Different hypergraph projections can be partitioned and used to create Partition objects that reflect the associated decomposition

abstract type AbstractPartition end

#A hypergraph community with a parent
struct PartitionLeaf <: AbstractPartition
    hypergraph::HyperGraph
    parent::Union{Nothing,AbstractPartition}
end

#PartitionParent descrives shared nodes and shared edges among its children
struct PartitionParent <: AbstractPartition
    sharednodes::Vector{HyperNode}          #shared nodes group into a master problem

    sharededges::Vector{HyperEdge}          #shared edges become link constraints

    parent::Union{Nothing,PartitionParent}
    children::Vector{AbstractPartition}
end
PartitionParent(sharednodes::Vector{HyperNode},sharededges::Vector{HyperEdge}) = PartitionParent(sharednodes,sharededges,nothing,Vector{AbstractPartition}())
PartitionParent(sharededges::Vector{HyperEdge}) = PartitionParent(Vector{HyperNode}(),sharededges,nothing,Vector{AbstractPartition}())

#A partition describes the entire partition structure.
mutable struct Partition
    subpartitions::Vector{PartitionLeaf}  #bottom level communities (i.e. Leaf Communities)
    parents::Vector{PartitionParent}          #tree structure describing recursive structure and shared nodes and edges
end
Partition() = Partition(Vector{PartitionLeaf}(),Vector{PartitionParent}())

# mutable struct OverlappingPartition
#     subpartitions::Vector{PartitionLeaf}
#     overlaps::Vector{Vector{HyperNode}}
# end

#Convert membership vector to lists of indices
function getpartitionlist(hypergraph::HyperGraph,membership_vector::Vector)
    unique_parts = unique(membership_vector)  #get unique membership entries
    unique_parts = sort(unique_parts)
    nparts = length(unique_parts)             #number of partitions

    partitions = OrderedDict{Int64,Vector{HyperNode}}((k,[]) for k in unique_parts)
    for (vertex,part) in enumerate(membership_vector)
        push!(partitions[part],getnode(hypergraph,vertex))
    end
    return collect(values(partitions))
end

function identifyhyperedges(hypergraph::HyperGraph,partitions::Vector{Vector{HyperNode}})
    nparts = length(partitions)

    #Create partition matrix
    I = []
    J = []
    for i = 1:nparts
       for hypernode in partitions[i]
           j = getindex(hypergraph,hypernode)
           push!(I,i)
           push!(J,j)
       end
    end

    V = Int.(ones(length(J)))
    G = sparse(I,J,V)  #Node partition matrix
    A = sparse(hypergraph)
    C = G*A  #Edge partitions

    #FIND THE SHARED EDGES, Get indices of shared edges
    sum_vector = sum(C,dims = 1)
    max_vector = maximum(C,dims = 1)
    cross_vector = sum_vector - max_vector
    indices = findall(cross_vector .!= 0)                   #nonzero indices of the cross vector.  These are edges that cross partitions.
    indices = [indices[i].I[2] for i = 1:length(indices)]   #convert to Integers

    shared_edges = HyperEdge[]
    for index in indices
        push!(shared_edges,gethyperedge(hypergraph,index))
    end

    #GET INDUCED PARTITION EDGES (I.E GET THE EDGES LOCAL TO EACH PARTITION)
    partition_edges = Vector[Vector{HyperEdge}() for _ = 1:nparts]
    for i = 1:nparts
        inds = findall(C[i,:] .!= 0)
        new_inds = filter(x -> !(x in indices), inds) #these are edge indices
        for new_ind in new_inds
            push!(partition_edges[i],gethyperedge(hypergraph,new_ind))
        end
    end

    return partition_edges,shared_edges
end

#Simple 2 level partition from a vector of integers
#Can be used for both a row-net HyperGraph or a clique-expansion Graph
function Partition(hypergraph::AbstractHyperGraph,node_membership_vector::Vector{Int64})
    hyperpartition = Partition()

    #convert membership vector to vector of vectors
    hypernode_vectors = getpartitionlist(hypergraph,node_membership_vector)

    #println("Identifying Shared and Induced Edges")
    induced_edge_partitions,shared_edges = identifyhyperedges(hypergraph,hypernode_vectors)

    #println("Creating Sub Hyper Graphs")
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
            add_sub_hyperedge!(hyper,hyperedge)  #NOTE HyperEdge Construction might be taking too long
        end

        push!(new_hypers,hyper)
    end

    #println("Creating Partition Parent")
    partition_parent = PartitionParent(shared_edges)
    partitions = Vector{SubgraphPartition}()
    for i = 1:length(new_hypers)
        push!(partitions,SubgraphPartition(new_hypers[i],partition_parent))
    end
    hyperpartition.partitions = partitions
    hyperpartition.parents = [partition_parent]

    return hyperpartition
end


#Create a ModelGraph Subgraph from a HyperGraph
function create_sub_modelgraph(modelgraph::ModelGraph,hypergraph::HyperGraph)
    submg = ModelGraph()
    submg.hypergraph = hypergraph

    for hypernode in getnodes(hypergraph)
        modelnode = getnode(modelgraph,hypernode)
        submg.modelnodes[hypernode] = modelnode
    end

    i = 1
    for hyperedge in getallhyperedges(hypergraph)
        linkedge = findlinkedge(modelgraph,hyperedge)  #could be in a subgraph
        submg.linkedges[hyperedge] = linkedge
        for linkconstraintref in linkedge.linkconstraints
            linkconstraint = LinkConstraint(linkconstraintref)
            submg.linkconstraints[i] = linkconstraint
            i += 1
        end
    end
    return submg
end



####################################
#Print Functions
####################################
function string(partition::Partition)
    """
    Partition:
    partitions: $(length(partition.partitions))
    """
end
print(io::IO, partition::Partition) = print(io, string(partition))
show(io::IO,partition::Partition) = print(io,partition)


#TODO
# function Partition(clique_graph::DualCliqueGraph,projection_map::ProjectionMap,membership_vector::Vector{Int64}) #NOTE: Could also be a Dual Clique Graph
#
#     hyperpartition = HyperPartition()
#
#     #figure out the hypergraph partition based on the graph partition
#
#     return hyperpartition
# end
#
# function Partition(bipartite_graph::BipartiteGraph,projection_map::ProjectionMap,membership_vector::Vector{Int64};selection = :shared_nodes)
#
#     hyperpartition = HyperPartition()
#
#     return hyperpartition
# end
#
# function Partition(dual_hyper_graph::AbstractHyperGraph,projection_map::ProjectionMap,membership_vector::Vector{Int64})
#
#     hyperpartition = HyperPartition()
#
#     return hyperpartition
# end
#

#NOTE (s)
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
