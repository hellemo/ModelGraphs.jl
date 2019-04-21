import CommunityDetection

"""
LightGraphs.label_propagation(graph::ModelGraph)

Return partitions corresponding to detected communities using the LightGraphs label propagation algorithm.
"""
function CommunityDetection.community_detection_nback(graph::ModelGraph,args...;kwargs...)
    return _getpartitions(graph,CommunityDetection.community_detection_nback,args...;kwargs...)
end

function CommunityDetection.community_detection_bethe(graph::ModelGraph,args...;kwargs...)
    return _getpartitions(graph,CommunityDetection.community_detection_bethe,args...;kwargs...)
end

function CommunityDetection.community_detection_louvain(graph::ModelGraph,args...;projection = NodeUnipartiteGraph,kwargs...)

    #NOTE: Return communities and shared entities.  These can be used to create an aggregated ModelGraph.
    return _getpartitions(graph,CommunityDetection.community_detection_louvain,projection = projection,args...;,kwargs...))
end
