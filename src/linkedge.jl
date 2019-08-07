##############################################################################
# ModelEdges
##############################################################################
struct LinkEdge <: AbstractLinkEdge
    hyperedge::HyperEdge                        #Reference to a HyperEdge
    linkconstraints::Vector{LinkConstraintRef}  #Link constraints this edge represents
end

function add_linkconstraint_edges!(graph::AbstractModelGraph,con_refs::Array{LinkConstraintRef}) #TODO make sure this always works
    for con_ref in con_refs
        add_link_edge!(graph,con_ref)
    end
end

function add_linkconstraint_edges!(graph::AbstractModelGraph,con_refs::JuMP.Containers.DenseAxisArray{AbstractGraphConstraintRef}) #TODO make sure this always works
    for con_ref in con_refs.data
        add_link_edge!(graph,con_ref)
    end
end

function add_link_edge!(graph::AbstractModelGraph,ref::LinkConstraintRef)

    linkconstraint = LinkConstraint(ref)   #Get the Linkconstraint object so we can inspect the nodes on it
    modelnodes = getnodes(linkconstraint)

    #Add hyper edge
    hypernodes = gethypernode.(modelnodes)
    hyper_edge = add_hyper_edge!(gethypergraph(graph),hypernodes...)

    #Map to LinkEdge
    #Either create new LinkEdge or look up existing one
    if haskey(graph.link_edges,hyper_edge)
        link_edge = graph.link_edges[hyper_edge]
    else
        link_edge = LinkEdge(hyper_edge)
        graph.linkedges[hyper_edge] = link_edge
    end
    push!(link_edge.linkconstraints,ref)
    graph.linkconstraint_linkedge_map[ref] = link_edge

    return link_edge
end


function string(edge::LinkEdge)
    "Linking edge w/ $(length(edge.linkconstraints)) Constraint(s)"
end
print(io::IO,edge::LinkingEdge) = print(io, string(edge))
show(io::IO,edge::LinkingEdge) = print(io,edge)



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
