#Create a ModelGraph from a given Hypergraph
function create_sub_modelgraph(modelgraph::ModelGraph,hypergraph::HyperGraph)
    submg = ModelGraph()
    submg.hypergraph = hypergraph

    for hypernode in getnodes(hypergraph)
        modelnode = getnode(modelgraph,hypernode)
        submg.modelnodes[hypernode] = modelnode
    end

    i = 1
    for hyperedge in getallhyperedges(hypergraph)
        linkedge = findlinkedge(modelgraph,hyperedge)  #this could be in a subgraph
        submg.linkedges[hyperedge] = linkedge
        for linkconstraintref in linkedge.linkconstraints
            linkconstraint = LinkConstraint(linkconstraintref)
            submg.linkconstraints[i] = linkconstraint
            i += 1
        end
    end
    return submg
end

function _create_worker_modelgraph(master::JuMP.Model,modelnodes::Vector{ModelNode},node_indices::Vector{Int64},n_nodes::Int64,n_linkeq_cons::Int64,n_linkineq_cons::Int64,linkeq_dict::OrderedDict,linkineq_dict::OrderedDict)
    graph = ModelGraph()
    graph.mastermodel = master

    #Add nodes to worker's graph.  Each worker should have the same number of nodes, but some will be empty.
    for i = 1:n_nodes
        add_node!(graph)
    end

    #Populate models for given nodes
    for (i,node) in enumerate(modelnodes)
        index = node_indices[i]  #need node index in highest level
        new_node = getnode(graph,index)
        set_model(new_node,getmodel(node))
        new_node.partial_linkconstraints = node.partial_linkconstraints
    end
    #We need the graph to have the partial constraints over graph nodes
    graph.linkconstraints = _add_link_terms(modelnodes)
    

    graph.obj_dict[:n_linkeq_cons] = n_linkeq_cons
    graph.obj_dict[:n_linkineq_cons] = n_linkineq_cons
    graph.obj_dict[:linkeq_dict] = linkeq_dict
    graph.obj_dict[:linkineq_dict] = linkineq_dict
    return graph
end

function _add_link_terms(modelnodes::Vector{ModelNode})
    linkconstraints = OrderedDict()
    for node in modelnodes
        partial_links = node.partial_linkconstraints
        for (idx,linkconstraint) in partial_links
            if !(haskey(linkconstraints,idx))   #create link constraint
                new_func = linkconstraint.func
                set = linkconstraint.set
                linkcon = LinkConstraint(new_func,set)
                linkconstraints[idx] = linkcon
            else #update linkconstraint
                newlinkcon = linkconstraints[idx]
                nodelinkcon = node.partial_linkconstraints[idx]
                newlinkcon = LinkConstraint(newlinkcon.func + nodelinkcon.func,newlinkcon.set)
                linkconstraints[idx] = newlinkcon
            end
        end
    end
    return linkconstraints
end
