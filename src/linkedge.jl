##############################################################################
# ModelEdges
##############################################################################
struct LinkEdge <: AbstractLinkEdge
    hyperedge::HyperEdge                        #Reference to a HyperEdge
    linkconstraints::Vector{AbstractLinkConstraintRef}  #Link constraints this edge represents
end
LinkEdge(hyperedge::HyperEdge) = LinkEdge(hyperedge,Vector{AbstractLinkConstraintRef}())


function add_link_edge!(graph::AbstractModelGraph,nodes::Vector{ModelNode})#ref::LinkConstraintRef)
    #Add hyper edge
    hypernodes = gethypernode.(graph,modelnodes)
    hyperedge = add_hyper_edge!(gethypergraph(graph),hypernodes...)

    #Map to LinkEdge
    #Either create new LinkEdge or look up existing one
    if haskey(graph.link_edges,hyperedge)
        link_edge = graph.link_edges[hyperedge]
    else
        link_edge = LinkEdge(hyper_edge)
        graph.linkedges[hyper_edge] = link_edge
    end

    push!(link_edge.linkconstraints,ref)

    graph.linkconstraint_linkedge_map[ref] = link_edge

    return link_edge
end


function string(edge::LinkEdge)
    "Link edge w/ $(length(edge.linkconstraints)) Constraint(s)"
end
print(io::IO,edge::LinkEdge) = print(io, string(edge))
show(io::IO,edge::LinkEdge) = print(io,edge)




# function add_linkconstraint_edges!(graph::AbstractModelGraph,con_refs::Array{LinkConstraintRef}) #TODO make sure this always works
#     for con_ref in con_refs
#         add_link_edge!(graph,con_ref)
#     end
# end
#
# function add_linkconstraint_edges!(graph::AbstractModelGraph,con_refs::JuMP.Containers.DenseAxisArray{AbstractGraphConstraintRef}) #TODO make sure this always works
#     for con_ref in con_refs.data
#         add_link_edge!(graph,con_ref)
#     end
# end


# Edge constructors
# LinkingEdge() = LinkingEdge(StructureGraphs.StructureEdge(),JuMP.ConstraintRef[])

# function add_edge!(graph::AbstractModelGraph)
#     add_hyper_edge!(graph.hypergraph)
#

#Add a hyperedge to graph using a linkconstraint reference
# function add_linkconstraint_edge!(graph::AbstractModelGraph,con_ref::LinkConstraintRef)
#
#     add_edge!(graph,con_ref)
#
#
# end
