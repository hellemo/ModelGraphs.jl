##############################################################################
# ModelEdges
##############################################################################
mutable struct LinkEdge <: AbstractLinkEdge
    vertices::Set{Int64}
    linkconstraints::Vector{AbstractLinkConstraintRef}  #Link constraints this edge represents
end
LinkEdge() = LinkEdge(Vector{Int64}(),Set{Int64}(),Vector{AbstractLinkConstraintRef}())

function add_link_edge!(graph::AbstractModelGraph,modelnodes::Vector{ModelNode}) #ref::LinkConstraintRef)
    node_indices = [graph.idx_map[node] for node in modelnodes]
    linkedge = LinkEdge(node_indices)
    n_links = length(graph.linkedges)
    idx = n_links + 1
    push!(graph.linkedges,linkedge)
    graph.linkedge_map[linkedge.vertices] = linkedge
    graph.edge_idx_map[linkedge] = idx
    return link_edge
end


function string(edge::LinkEdge)
    "Link edge w/ $(length(edge.linkconstraints)) Constraint(s)"
end
print(io::IO,edge::LinkEdge) = print(io, string(edge))
show(io::IO,edge::LinkEdge) = print(io,edge)
