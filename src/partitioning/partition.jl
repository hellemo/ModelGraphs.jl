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
function partition(graph::ModelGraph,projection = NodeUnipartiteGraph,partition_func::Function = Metis.partition,args...;kwargs...)

    projected_graph,projection_map = projection(graph)          #Reference map maps projected graph entities to model_graph entities

    partitions,local_entities,shared_entities = partition(projected_graph,partition_func,args...;kwargs...)

    #now map partition data back to original graph
    return_partitions = _map_partitions(partitions,projection_map)

    return_shared_entities = unique(_map_entities(shared_entities,projection_map))
    #NOTE: Need to keep vector the same size. #If there are duplicate entries across partitions, then they must also show up in shared
    return_partition_entities = [unique(_map_entities(local_entitiy,projection_map)) for local_entitiy in local_entities]

    #Make sure no partition entities are in shared entities.  It's possible that a local entity maps to a linkconstraint that is actually shared.
    for return_part in return_partition_entities
        filter!(e ->  !(e in return_shared_entities),return_part)
    end

    return PartitionData(return_partitions,return_partition_entities,return_shared_entities)
end

#return partition data
function partition(ugraph::NodeUnipartiteGraph,partition_func::Function,args...;kwargs...)
    #lg = getlightgraph(ugraph)
    #membership_vector = partition_func(lg,args...;kwargs...)
    membership_vector = partition_func(ugraph,args...;kwargs...)

    partitions = _getpartitions(ugraph,membership_vector)

    local_entities,shared_entities = _identifyentities(ugraph,partitions)  #Should be vector of edges in the graph

    return partitions,local_entities,shared_entities

    # return_partitions = _map_partitions(partitions,projection_map)
    # return_shared_entities = unique(_map_entities(shared_entities,projection_map))
    #
    # #NOTE: Need to keep vector the same size. #If there are duplicate entries across partitions, then they must also show up in shared
    # return_partition_entities = [unique(_map_entities(local_entitiy,projection_map)) for local_entitiy in local_entities]
    #
    # #Make sure no partition entities are in shared entities.  It's possible that a local entity maps to a linkconstraint that is actually shared.
    # for return_part in return_partition_entities
    #     filter!(e ->  !(e in return_shared_entities),return_part)
    # end
    #
    # return PartitionData(return_partitions,return_partition_entities,return_shared_entities)#typeof(ugraph))
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

    return partitions
end

#Get shared linkconstraints between partitions
#TODO Take a LinearAlgebra approach for this (see some reference papers on how to do this)
function _identifyentities(graph::NodeUnipartiteGraph,partitions::Vector{Vector{Int64}})
    println("Identifying Shared entities")
    nparts = length(partitions)

    #FIND THE SHARED EDGES
    #Setup sparse partition matrix
    I = []
    J = []
    for i = 1:nparts
        for j in partitions[i]
            push!(I,i)
            push!(J,j)
        end
    end
    V = ones(length(J))
    G = sparse(I,J,V)  #Node partition matrix
    lg = getlightgraph(graph)
    A = LightGraphs.incidence_matrix(lg)
    C = G*A  #Edge partitions

    #Get indices of shared edges
    sum_vector = sum(C,dims = 1)
    max_vector = maximum(C,dims = 1)
    cross_vector = sum_vector - max_vector
    indices = findall(cross_vector .!= 0)
    indices = [indices[i].I[2] for i = 1:length(indices)]  #some data unpacking to get the indices

    cross_matrix = A[:,indices]
    n_cross_edges = size(cross_matrix)[2]
    shared_edges = LightGraphs.Edge[]
    #TODO: Need to check edges in both directions.  This change would be fundamental to my graph interface
    for i = 1:n_cross_edges
        node_indices = cross_matrix[:,i].nzind
        push!(shared_edges,LightGraphs.Edge(node_indices[1],node_indices[2]))
    end

    #GET THE EDGES LOCAL TO EACH PARTITION
     partition_edges = Vector[Vector{LightGraphs.Edge}() for _ = 1:nparts]
     for i = 1:nparts
         inds = findall(C[i,:] .!= 0)
         new_inds = filter(x -> !(x in indices), inds)
         local_matrix = A[:,new_inds]
         for j = 1:length(new_inds)
             node_indices = local_matrix[:,j].nzind
             push!(partition_edges[i],LightGraphs.Edge(node_indices[1],node_indices[2]))
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
function _map_entities(entities::Vector{Any},projection_map::AlgebraicGraphs.ProjectionMap)
    mapped_entities = []
    for entity in entities
        mapped_indices = projection_map[entity]   #this is a vector
        append!(mapped_entities,mapped_indices)
        #push!(mapped_entities,mapped_index)
    end
    return mapped_entities
end

function _map_entities(edge_list::Vector{LightGraphs.SimpleGraphs.SimpleEdge},projection_map::AlgebraicGraphs.ProjectionMap)
    mapped_entities = []
    for edge in edge_list
        try
            # mapped_index = projection_map[edge]
            # push!(mapped_entities,mapped_index)
            mapped_indices = projection_map[edge]
            append!(mapped_entities,mapped_indices)
        catch KeyError
            rev_edge = LightGraphs.Edge(edge.dst,edge.src)
            # mapped_index = projection_map[rev_edge]
            # push!(mapped_entities,mapped_index)
            mapped_indices = projection_map[rev_edge]
            append!(mapped_entities,mapped_indices)
        end
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

#function _identifyentities2(graph::NodeUnipartiteGraph,partitions::Vector{Vector{Int64}})
    #println("Identifying Shared entities")
    # n_partitions = length(partitions)
    # partition_edges = Vector[Vector() for _ = 1:n_partitions]
    # shared_edges = []
    # checked_edges = []
    # for i = 1:n_partitions
    #     partition = partitions[i]
    #     #NOTE: We effectively need to iterate the edges and check for supporting nodes that are outside the partition
    #     for node_index in partition
    #         node = getnode(graph,node_index)
    #         edges = getsupportingedges(graph,node)
    #         for edge in edges
    #             if !(edge in checked_edges)
    #                 edge_nodes = getsupportingnodes(graph,edge)
    #                 node_indices = [getindex(graph,node) for node in edge_nodes]
    #                 if !(all(node -> node in partition,node_indices))
    #                     push!(shared_edges,getindex(graph,edge))
    #                 else
    #                     push!(partition_edges[i],getindex(graph,edge))
    #                 end
    #             end
    #             push!(checked_edges,edge)
    #         end
    #     end
    # end
    #return partition_edges,shared_edges
#end
