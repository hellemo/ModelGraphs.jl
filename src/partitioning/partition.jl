const PartitionEntity = Union{ModelNode,LinkingEdge,Pair{ModelNode,LinkingEdge}}
#Partition data pertaining to partitions of nodes, edges, or pairs depending on the projection used to create it
# struct PartitionData
#     partitions::Vector{Vector{Union{ModelNode,LinkingEdge}}}  #partitions of nodes or edges or both
#     partition_entities::Vector{Vector}
#     shared_entities::Vector{PartitionEntity}
#     partition_type::Function  #Nodes, Edges, or Pairs
# end
#TODO Typing
struct PartitionData
    partitions::Vector{Vector}  #partitions of nodes or edges or both
    partition_entities::Vector{Vector}
    shared_entities::Vector{}
    #partition_type::Function  #Nodes, Edges, or Pairs
end

#TODO
function get_shared_entities(partition_data::PartitionData)
end


#Helper function to convert vector of node membership in partitions into vectors of node indices
function partition(graph::ModelGraph,partition_func::Function,projection = NodeUnipartiteGraph,args...;kwargs...)
    projected_graph,projection_map = projection(graph)          #Reference map maps projected graph entities to model_graph entities
    partition_data = partition(projected_graph,partition_func,projection_map,args...;kwargs...)
    return partition_data
end

#return partition data
function partition(ugraph::NodeUnipartiteGraph,partition_func::Function,projection_map::ProjectionMap,args...;kwargs...)
    lg = getlightgraph(ugraph)
    membership_vector = partition_func(lg,args...;kwargs...)

    partitions = _getpartitions(ugraph,membership_vector)
    local_entities,shared_entities = _identifyentities(ugraph,partitions)  #Should be vector of edges in the graph

    return_partitions = _map_partitions(partitions,projection_map)
    return_shared_entities = unique(_map_entities(shared_entities,projection_map))
    return_partition_entities = unique(map(x -> _map_entities(x,projection_map),local_entities))

    return PartitionData(return_partitions,return_partition_entities,return_shared_entities)#typeof(ugraph))
end

#Convert membership vector to lists of indices
function _getpartitions(graph::NodeUnipartiteGraph,membership_vector::Vector)
    unique_parts = unique(membership_vector)            #get unique membership entries
    nparts = length(unique_parts)           #number of partitions
    partition_dict = Dict{Int64,Vector{Int64}}((k,[]) for k in unique_parts)

    for node in getnodes(graph)
        index = getindex(graph,node)  #index in the UnipartiteGraph
        part = membership_vector[index]      #parition id of this node index
        push!(partition_dict[part],index)
    end

    partitions = collect(values(partition_dict))  #returns partitions of original model graph nodes

    # NOTE: A way to get a sparse group identity matrix
    # make sure graph has no self loops
    # I = []
    # J = []
    # for i = 1:nparts
    #     for j in partitions[i]
    #         push!(I,i)
    #         push!(J,j)
    #     end
    # end
    # V = ones(length(J))
    # G = sparse(I,J,V)

    return partitions
end

#Get shared linkconstraints between partitions
#TODO Take a LinearAlgebra approach for this (see some reference papers on how to do this)
#NOTE: This function is too slow.  Look at other approaches to find shared entities
#NOTE: Need to fix hypergraph implementation so I don't get self loops in the light graph
function _identifyentities(graph::NodeUnipartiteGraph,partitions::Vector{Vector{Int64}})
    println("Identifying Shared entities")
    n_partitions = length(partitions)
    partition_edges = Vector[Vector() for _ = 1:n_partitions]
    shared_edges = []
    checked_edges = []
    for i = 1:n_partitions
        partition = partitions[i]
        #NOTE: We effectively need to iterate the edges and check for supporting nodes that are outside the partition
        for node_index in partition
            node = getnode(graph,node_index)
            edges = getsupportingedges(graph,node)
            for edge in edges
                if !(edge in checked_edges)
                    edge_nodes = getsupportingnodes(graph,edge)
                    node_indices = [getindex(graph,node) for node in edge_nodes]
                    if !(all(node -> node in partition,node_indices))
                        push!(shared_edges,getindex(graph,edge))
                    else
                        push!(partition_edges[i],getindex(graph,edge))
                    end
                end
                push!(checked_edges,edge)
            end
        end
    end
    return partition_edges,shared_edges
end

#TODO Typing output
function _map_partitions(partitions::Vector{Vector{Int64}},projection_map::AlgebraicGraphs.ProjectionMap)
    n_partitions = length(partitions)
    mapped_partitions = Vector[Vector() for _ = 1:n_partitions]

    for i = 1:n_partitions
        partition = partitions[i]
        for index in partition
            mapped_index = projection_map[index]
            push!(mapped_partitions[i],mapped_index)
        end
    end
    return mapped_partitions
end


#TODO Typing input and output
function _map_entities(entities::Vector,projection_map::AlgebraicGraphs.ProjectionMap)
    mapped_entities = []
    for entity in entities
        mapped_index = projection_map[entity]
        push!(mapped_entities,mapped_index)
    end
    return mapped_entities
end

# function _map_entities(entities::Vector{Vector{LightGraphs.AbstractEdge}},projection_map::ProjectionMap)
#     for vector in entities
#
#     mapped_entities = []
#     for entity in entities
#         mapped_index = projection_map[entity]
#         push!(mapped_entities,mapped_index)
#     end
# end


# #NOTE: Trying to update how this works
# function _get_local_and_cross_links(model_graph::ModelGraph,nodes::Vector{ModelNode})
#     local_link_constraints = []  #all nodes in link constraint are in the set of nodes given
#     cross_link_constraints = []  #at least one node in a link constraint is not in the node subset (must be connected to the rest of the graph)
#     #Check each node's link constraints.  Add it to one of the above lists.
#     checked_links = []
#     for node in nodes
#         for (graph,links) in getlinkconstraints(node)  #NOTE:  Need to store linkconstraint references on nodes
#             for link in links
#                 if !(link in checked_links)  #If it's a new link constraint
#                     vars = link.terms.vars
#                     var_nodes = map(getnode,vars)
#                     if all(node -> node in nodes,var_nodes)
#                         push!(local_link_constraints,link)
#                     else
#                         push!(cross_link_constraints,link)
#                     end
#                     push!(checked_links,link)
#                 end
#             end
#         end
#     end
#     return local_link_constraints , cross_link_constraints
# end
