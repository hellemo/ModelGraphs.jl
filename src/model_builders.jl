#Create a ModelGraph Subgraph from a HyperGraph
function create_sub_modelgraph(modelgraph::ModelGraph,hypergraph::HyperGraph)
    submg = ModelGraph()
    submg.hypergraph = hypergraph

    for hypernode in getnodes(hypergraph)
        modelnode = getnode(modelgraph,hypernode)
        submg.modelnodes[hypernode] = modelnode
    end

    i = 1
    for hyperedge in getallhyperedges(hypergraph)
        linkedge = findlinkedge(modelgraph,hyperedge)  #could be in a subgraph
        submg.linkedges[hyperedge] = linkedge
        for linkconstraintref in linkedge.linkconstraints
            linkconstraint = LinkConstraint(linkconstraintref)
            submg.linkconstraints[i] = linkconstraint
            i += 1
        end
    end
    return submg
end


function _create_worker_modelgraph(master::JuMP.Model,modelnodes::Vector{ModelNode},n_nodes::Int64,n_linkeq_cons::Int64,n_linkineq_cons::Int64)
    graph = ModelGraph()
    graph.master = master

    #Add nodes to worker's graph.  Each worker should have the same number of nodes, but some will be empty.
    for i = 1:n_nodes
        add_node!(graph)
    end

    #Populate models for given nodes
    for node in modelnodes
        index = node.ext[:index]  #need node index in highest level
        new_node = getnode(graph,index)
        setmodel(new_node,getmodel(node))
        new_node.partial_linkconstraints = node.partial_linkconstraints
        #copy_partial_linkconstraints!(node,new_node)  #create copies of partial linkconstraints
    end

    #We need the graph to have the partial constraints over graph nodes
    graph.linkconstraints = _add_link_terms(modelnodes)
    graph.obj_dict[:n_linkeq_cons] = n_linkeq_cons
    graph.obj_dict[:n_linkineq_cons] = n_linkineq_cons 
    return graph
end

function _add_link_terms(modelnodes::Vector{ModelNode})
    linkconstraints = Dict()
    for node in modelnodes
        partial_links = node.partial_linkconstraints
        for (idx,linkconstraint) in partial_links
            if !(haskey(linkconstraints,idx))   #create link constraint
                new_func = linkconstraint.func
                set = linkconstraint.set
                linkcon = LinkConstraint(new_func,set)
                linkconstraints[idx] = linkcon
            else #update linkconstraint
                newlinkcon = node.partial_linkconstraints[index]
                JuMP.add_to_expression!(newlinkcon.func,linkconstraint.func)
            end
        end
    end
    return linkconstraints
end
