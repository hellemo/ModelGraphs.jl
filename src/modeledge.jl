##############################################################################
# Edges
##############################################################################
struct LinkingEdge <: AbstractLinkingEdge
    linkedge::StructureGraphs.StructureEdge
    linkconstraints::Vector{ConstraintRef}  #Link constraints this edge represents
end
#Edge constructors
LinkingEdge() = LinkingEdge(StructureGraphs.StructureEdge(),JuMP.ConstraintRef[])
create_edge(graph::ModelGraph) = LinkingEdge()

#Add a hyperedge to graph using a linkconstraint reference
function add_edge!(graph::AbstractModelGraph,ref::GraphConstraintRef)

    con = LinkConstraint(ref)   #Get the Linkconstraint object so we can inspect the nodes on it
    nodes = getnodes(con)

    #node_indices = con.node_indices

    edge = add_edge!(graph,nodes...)
    push!(edge.linkconstraints,ref)

    #NOTE: Is storing this information necessary?.  We can look at a node's incident edges and determine the linkconstraints.

    # for node in nodes
    #     push!(node.linkconrefs[graph],ref)
    # end

    # for node in nodes
    #     if !haskey(node.linkconstraints,graph)
    #         node.linkconstraints[graph] = [ref]
    #     else
    #         push!(node.linkconrefs[graph],ref)
    #     end
    # end

    return edge
end

# TODO  Think of a good way to update links when swapping out models.  Might need to store variable names in NodeLinkData
# function _updatelinks(m,::AbstractModel,nodeoredge::NodeOrEdge)
#     link_cons = getlinkconstraints(nodeoredge)
#     #find variables
# end
