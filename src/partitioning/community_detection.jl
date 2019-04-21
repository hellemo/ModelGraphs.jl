"""
LightGraphs.label_propagation(graph::ModelGraph)

Return partitions corresponding to detected communities using the LightGraphs label propagation algorithm.
"""
function LightGraphs.label_propagation(graph::ModelGraph,projection::Function = NodeUnipartiteGraph)
    return _getpartitions(graph,LightGraphs.label_propagation,projection)
end
