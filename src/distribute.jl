#Distribute a modelgraph among workers.  Each worker should have the same master model.  Each worker will be allocated some of the nodes in the modelgraph
function distribute(mg::ModelGraph,to_workers::Vector{Int64};remote_name = :graph)
    #NOTE: Linkconstraints keep their indices in new graphs, NOTE: Link constraint row index needs to match on each worker
    #NOTE: Does not yet support subgraphs.  Aggregate first
    #Create remote channel to store the nodes we want to send

    #IDEA: Create a channel from the master process to each worker?

    channel_nodes = RemoteChannel(1)    #we will allocate and send nodes to workers
    channel_indices = RemoteChannel(1)
    channel_master = RemoteChannel(1)   #we will send the master problem to each worker

    n_nodes = getnumnodes(mg)
    n_workers = length(to_workers)
    nodes_per_worker = Int64(floor(n_nodes/n_workers))
    nodes = all_nodes(mg)
    node_indices = [getindex(mg,node) for node in nodes]

    #link_data = ModelGraphs.get_link_constraint_data(mg)
    n_linkeq_cons = length(mg.linkeqconstraints)
    n_linkineq_cons = length(mg.linkineqconstraints)
    #n_linkineq_cons = length(link_data.linear_le_constraints) + length(link_data.linear_ge_constraints) + length(link_data.linear_interval_constraints)
    n_link_cons = n_linkeq_cons + n_linkineq_cons

    #Allocate modelnodes onto provided workers
    allocations = []
    node_indices = []
    j = 1
    while  j <= n_nodes
        if j + nodes_per_worker > n_nodes
            push!(allocations,nodes[j:end])
            push!(node_indices,[getindex(mg,node) for node in nodes[j:end]])
        else
            push!(allocations,nodes[j:j+nodes_per_worker - 1])
            push!(node_indices,[getindex(mg,node) for node in nodes[j:j+nodes_per_worker - 1]])
        end
        j += nodes_per_worker
    end
    master = getmasternode(mg)
    put!(channel_master, [master])  #put master model (node) into channel

    println("Distributing graph among workers: $to_workers")
    remote_references = []
    #Fill channel with sets of nodes to send
    @sync begin
        for (i,worker) in enumerate(to_workers)
            @spawnat(1, put!(channel_nodes, allocations[i]))
            @spawnat(1, put!(channel_indices, node_indices[i]))
            ref1 = @spawnat worker begin
                Core.eval(Main, Expr(:(=), :master, fetch(channel_master)[1]))
                Core.eval(Main, Expr(:(=), :nodes, take!(channel_nodes)))
                Core.eval(Main, Expr(:(=), :node_indices, take!(channel_indices)))
            end
            wait(ref1)
            ref2 = @spawnat worker Core.eval(Main, Expr(:(=), remote_name, ModelGraphs._create_worker_modelgraph(getfield(Main,:master),getfield(Main,:nodes),getfield(Main,:node_indices),
            n_nodes,n_linkeq_cons,n_linkineq_cons,mg.linkeqconstraints,mg.linkineqconstraints)))
            push!(remote_references,ref2)
        end
        return remote_references
    end
end

function _create_worker_modelgraph(master::ModelNode,modelnodes::Vector{ModelNode},node_indices::Vector{Int64},n_nodes::Int64,n_linkeq_cons::Int64,n_linkineq_cons::Int64,linkeq_dict::OrderedDict,linkineq_dict::OrderedDict)
    graph = ModelGraph()
    graph.node_idx_map = Dict{ModelNode,Int64}()
    graph.masternode = master
    graph.node_idx_map[master] = 0

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
    #graph.linkconstraints = _add_link_terms(modelnodes)
    graph.linkeqconstraints = _add_linkeq_terms(modelnodes)
    graph.linkineqconstraints = _add_linkineq_terms(modelnodes)

    #Tell the worker how many linkconstraints the graph actually has
    graph.obj_dict[:n_linkeq_cons] = n_linkeq_cons
    graph.obj_dict[:n_linkineq_cons] = n_linkineq_cons
    #graph.obj_dict[:linkeq_dict] = linkeq_dict
    #graph.obj_dict[:linkineq_dict] = linkineq_dict
    return graph
end

function _add_linkeq_terms(modelnodes::Vector{ModelNode})
    linkeqconstraints = OrderedDict()
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
