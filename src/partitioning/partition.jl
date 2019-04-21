const PartitionEntity = Union{ModelNode,LinkingEdge,Pair{ModelNode,LinkingEdge}}
#Partition data pertaining to partitions of nodes, edges, or pairs depending on the projection used to create it
struct PartitionData
    partitions::Vector{Vector{Union{ModelNode,LinkingEdge}}}  #partitions of nodes or edges or both
    shared_entities::Vector{PartitionEntity}
    partition_type::Function  #Nodes, Edges, or Pairs
end

#Helper function to convert vector of node membership in partitions into vectors of node indices
function partition(graph::ModelGraph,partition_func::Function,projection::Function = NodeUnipartiteGraph,args...;kwargs...)
    projected_graph,projection_map = projection(graph)          #Reference map maps projected graph entities to model_graph entities
    partition_data = partition(projected_graph,partition_func,projection_map,args...;kwargs...)
    return partition_data
end

#return partition data
function partition(ugraph::NodeUnipartiteGraph,partition_func::Function,projection_map::ProjectionMap,args...;kwargs...)
    lg = getlightgraph(ugraph)
    membership_vector = partition_func(lg,args...;kwargs...)

    partitions = _getpartitions(ugraph,projection_map,membership_vector)
    shared_entities = _getsharedentities(ugraph,partitions)  #Should be vector of edges in the graph

    return_partitions = _map_partitions(partitions,projection_map)
    return_shared_entities = _map_shared_entities(shared_entities,projection_map)

    return PartitionData(return_partitions,return_shared_entities)
end

function _getpartitions(graph::NodeUnipartiteGraph,membership_vector::Vector)
    unique_parts = unique(parts)            #get unique membership entries
    nparts = length(unique_parts)           #number of partitions
    partition_dict = Dict{Int64,Vector{Int64}}((k,[]) for k in unique_parts)

    for node in getnodes(graph)
        index = getindex(graph,node)  #index in the UnipartiteGraph
        part = membership_vector[index]      #parition id of this node index

        #model_node_index = projection_map[index]       #index of the node
        push!(partition_dict[part],index)
    end

    partitions = collect(values(partition_dict))  #returns partitions of original model graph nodes
    return partitions
end

#Get shared linkconstraints between partitions
#TODO Take a LinearAlgebra approach for this (see some reference papers on how to do this)
function _getsharedentities(graph::NodeUnipartiteGraph,partitions::Vector{Vector{Int64}})
    shared_edges = []
    checked_edges = []
    for partition in partitions
        for node_index in partition
            node = getnode(graph,node_index)
            edges = getincidentedges(graph,node)
            for edge in edges
                if !(edge in checked_edges)
                    edge_nodes = getsupportingnodes(graph,edge)
                    if !(all(node -> node in partition,edge_nodes))
                        push!(shared_edges,edge)
                    end
                end
                push!(checked_edges,edge)
            end
        end
    end
end







        #     for (graph,links) in getlinkconstraints(node)  #NOTE:  Need to store linkconstraint references on nodes
        #         for link in links
        #             if !(link in checked_links)  #If it's a new link constraint
        #                 vars = link.terms.vars
        #                 var_nodes = map(getnode,vars)
        #                 if all(node -> node in nodes,var_nodes)
        #                     push!(local_link_constraints,link)
        #                 else
        #                     push!(cross_link_constraints,link)
        #                 end
        #                 push!(checked_links,link)
        #             end
        #         end
        #     end
        # end


end


#NOTE: Trying to update how this works
function _get_local_and_cross_links(model_graph::ModelGraph,nodes::Vector{ModelNode})
    local_link_constraints = []  #all nodes in link constraint are in the set of nodes given
    cross_link_constraints = []  #at least one node in a link constraint is not in the node subset (must be connected to the rest of the graph)
    #Check each node's link constraints.  Add it to one of the above lists.
    checked_links = []
    for node in nodes
        for (graph,links) in getlinkconstraints(node)  #NOTE:  Need to store linkconstraint references on nodes
            for link in links
                if !(link in checked_links)  #If it's a new link constraint
                    vars = link.terms.vars
                    var_nodes = map(getnode,vars)
                    if all(node -> node in nodes,var_nodes)
                        push!(local_link_constraints,link)
                    else
                        push!(cross_link_constraints,link)
                    end
                    push!(checked_links,link)
                end
            end
        end
    end
    return local_link_constraints , cross_link_constraints
end
