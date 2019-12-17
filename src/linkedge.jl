##############################################################################
# ModelEdges
##############################################################################
struct LinkEdge <: AbstractLinkEdge
    hyperedge::HyperEdge                        #Reference to a HyperEdge
    linkconstraints::Vector{AbstractLinkConstraintRef}  #Link constraints this edge represents
end
LinkEdge(hyperedge::HyperEdge) = LinkEdge(hyperedge,Vector{AbstractLinkConstraintRef}())

function add_link_edge!(graph::AbstractModelGraph,modelnodes::Vector{ModelNode})#ref::LinkConstraintRef)
    #Add hyper edge
    hypernodes = gethypernode.(modelnodes)
    hyperedge = add_hyperedge!(gethypergraph(graph),hypernodes...)

    #Map to LinkEdge
    #Either create new LinkEdge or look up existing one
    if haskey(graph.linkedges,hyperedge)
        link_edge = graph.linkedges[hyperedge]
    else
        link_edge = LinkEdge(hyperedge)
        graph.linkedges[hyperedge] = link_edge
    end
    return link_edge
end


function string(edge::LinkEdge)
    "Link edge w/ $(length(edge.linkconstraints)) Constraint(s)"
end
print(io::IO,edge::LinkEdge) = print(io, string(edge))
show(io::IO,edge::LinkEdge) = print(io,edge)
