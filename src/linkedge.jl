##############################################################################
# ModelEdges
##############################################################################
struct LinkingEdge <: AbstractLinkingEdge
    hyperedge::StructuredGraphs.HyperEdge
    linkconstraints::Vector{AbstractGraphConstraintRef}  #Link constraints this edge represents
end
#Edge constructors
# LinkingEdge() = LinkingEdge(StructureGraphs.StructureEdge(),JuMP.ConstraintRef[])
# StructureGraphs.create_edge(graph::AbstractModelGraph) = LinkingEdge()
# StructureGraphs.getstructureedge(edge::LinkingEdge) = edge.structureedge

#Add a hyperedge to graph using a linkconstraint reference
function addlinkconstraintedge!(graph::AbstractModelGraph,con_ref::AbstractGraphConstraintRef)
    add_edge!(graph,con_ref)
end

function addlinkconstraintedges!(graph::AbstractModelGraph,con_refs::Array{AbstractGraphConstraintRef}) #TODO make sure this always works
    for con_ref in con_refs
        add_edge!(graph,con_ref)
    end
end

function addlinkconstraintedges!(graph::AbstractModelGraph,con_refs::JuMP.Containers.DenseAxisArray{AbstractGraphConstraintRef}) #TODO make sure this always works
    for con_ref in con_refs.data
        add_edge!(graph,con_ref)
    end
end

function add_hyper_edge!(graph::AbstractModelGraph,ref::AbstractGraphConstraintRef)

    con = LinkConstraint(ref)   #Get the Linkconstraint object so we can inspect the nodes on it
    nodes = getnodes(con)

    #node_indices = con.node_indices

    add_hyper_edge!(graph.hypergraph,node_indices)
    #edge = StructureGraphs.add_edge!(graph,nodes...)   #Add a new hyperedge given the nodes

    push!(edge.linkconstraints,ref)

    #NOTE: Consider storing linkconstraint references on nodes
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


function string(edge::LinkingEdge)
    "Linking edge w/ $(length(edge.linkconstraints)) Constraint(s)"
end
print(io::IO,edge::LinkingEdge) = print(io, string(edge))
show(io::IO,edge::LinkingEdge) = print(io,edge)
