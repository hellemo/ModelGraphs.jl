##############################################################################
# LinkEdges
# LinkEdges describe connections between model nodes
##############################################################################
mutable struct LinkEdge <: AbstractLinkEdge
    #vertices::Set{Int64}
    nodes::OrderedSet{ModelNode}
    linkconstraints::Vector{AbstractLinkConstraintRef}  #Link constraints this edge represents
end
LinkEdge() = LinkEdge(OrderedSet{ModelNode}(),Vector{AbstractLinkConstraintRef}())
LinkEdge(nodes::Vector{ModelNode}) = LinkEdge(OrderedSet(nodes),Vector{AbstractLinkConstraintRef}())
#LinkEdge() = LinkEdge(Set{Int64}(),Vector{AbstractLinkConstraintRef}())
#LinkEdge(indices::Set{Int64}) = LinkEdge(indices,Vector{AbstractLinkConstraintRef}())

function string(edge::LinkEdge)
    "Link edge w/ $(length(edge.linkconstraints)) Constraint(s)"
end
print(io::IO,edge::LinkEdge) = print(io, string(edge))
show(io::IO,edge::LinkEdge) = print(io,edge)


# function add_link_edge!(graph::AbstractModelGraph,modelnodes::Vector{ModelNode})
#     node_indices = Set([graph.node_idx_map[node] for node in modelnodes])
#     linkedge = LinkEdge(node_indices)
#     n_links = length(graph.linkedges)
#     idx = n_links + 1
#     push!(graph.linkedges,linkedge)
#     graph.linkedge_map[linkedge.vertices] = linkedge
#     graph.edge_idx_map[linkedge] = idx
#     return linkedge
# end
