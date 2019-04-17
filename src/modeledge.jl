##############################################################################
# ModelEdges
##############################################################################
struct LinkingEdge <: AbstractLinkingEdge
    structureedge::StructureGraphs.StructureEdge
    linkconstraints::Vector{GraphConstraintRef}  #Link constraints this edge represents
end
#Edge constructors
LinkingEdge() = LinkingEdge(StructureGraphs.StructureEdge(),JuMP.ConstraintRef[])
StructureGraphs.create_edge(graph::ModelGraph) = LinkingEdge()
StructureGraphs.getstructureedge(edge::LinkingEdge) = edge.structureedge

#Add a hyperedge to graph using a linkconstraint reference
function addlinkedges!(graph::AbstractModelGraph,con_ref::GraphConstraintRef)
    add_edge!(graph,con_ref)
end

function addlinkedges!(graph::AbstractModelGraph,con_refs::Array{GraphConstraintRef}) #TODO make sure this always works
    for con_ref in con_refs
        add_edge!(graph,con_ref)
    end
end

function addlinkedges!(graph::AbstractModelGraph,con_refs::JuMP.Containers.DenseAxisArray{GraphConstraintRef}) #TODO make sure this always works
    for con_ref in con_refs.data
        add_edge!(graph,con_ref)
    end
end

function StructureGraphs.add_edge!(graph::AbstractModelGraph,ref::GraphConstraintRef)

    con = LinkConstraint(ref)   #Get the Linkconstraint object so we can inspect the nodes on it
    nodes = getnodes(con)

    #node_indices = con.node_indices

    edge = StructureGraphs.add_edge!(graph,nodes...)
    push!(edge.linkconstraints,ref)




    #STORE LINKCONSTRAINT REFERENCES ON NODES
    #NOTE: Is storing this information necessary?.  We can look at a node's incident edges and determine the linkconstraints.
    # Depends how often we need to look up this information.  Could cache the references onto the nodes.
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
