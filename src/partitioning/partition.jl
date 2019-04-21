const PartitionEntity = Union{ModelNode,LinkingEdge,Pair{ModelNode,LinkingEdge}}

#Partition data pertaining to partitions of nodes, edges, or pairs depending on the projection used to create it
struct PartitionData
    partitions::Vector{Vector{Union{ModelNode,LinkingEdge}}}  #partitions of nodes or edges or both
    shared_entities::Vector{PartitionEntity}
    partition_type::Function  #Nodes, Edges, or Pairs
end

#Helper function to convert vector of node membership in partitions into vectors of node indices
function partition(graph::ModelGraph,partition_func::Function,projection::Function = NodeUnipartiteGraph,args...;kwargs...)
    projected_graph,reference_map = projection(graph) #Reference map maps projected graph entities to model_graph entities

    partition_data = partition_func(projected_graph,args...;kwargs...)

    #shared_entities = get_shared_entities(projected_graph,partitions)  #pass reference map here?

    return PartitionData(partitions,shared_entities,projection)
end

#return partition data
function partition(ugraph::NodeUnipartiteGraph,partition_func::Function,args...;kwargs...)
    lg = getlightgraph(ugraph)

    membership = partition_func(lg,args...;kwargs...)

    partitions = _getpartitions(membership)
    shared_entities = _getsharedentities(ugraph,partitions)

    return PartitionData(partitions,shared_entities)
end

function _getpartitions(membership::Vector,ref_map::ReferenceMap)
    unique_parts = unique(parts)
    nparts = length(unique_parts)
    partition_dict = Dict{Int64,Vector{Int64}}((k,[]) for k in unique_parts)

    for modelnode in getnodes(graph)
        index = getindex(graph,modelnode)
        part = parts[index]
        push!(partition_dict[part],index)
    end

    partitions = collect(values(partition_dict))
    return partitions
end
