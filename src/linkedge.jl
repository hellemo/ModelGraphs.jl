##############################################################################
# LinkEdges
# LinkEdges describe connections between model nodes
##############################################################################
mutable struct LinkEdge <: AbstractLinkEdge
    nodes::OrderedSet{ModelNode}
    linkconstraints::Vector{AbstractLinkConstraintRef}  #Link constraints this edge represents
    dual_values::Dict{AbstractLinkConstraintRef,Float64}
end
LinkEdge() = LinkEdge(OrderedSet{ModelNode}(),Vector{AbstractLinkConstraintRef}(),Dict{AbstractLinkConstraintRef,Float64}())
LinkEdge(nodes::Vector{ModelNode}) = LinkEdge(OrderedSet(nodes),Vector{AbstractLinkConstraintRef}(),Dict{AbstractLinkConstraintRef,Float64}())

function string(edge::LinkEdge)
    "Link edge w/ $(length(edge.linkconstraints)) Constraint(s)"
end
print(io::IO,edge::LinkEdge) = print(io, string(edge))
show(io::IO,edge::LinkEdge) = print(io,edge)
