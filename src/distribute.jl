#Distribute a modelgraph among workers.  Each worker should have the same master model.  Each worker will be allocated some of the nodes in the modelgraph
function ModelGraphs.distribute(mg::ModelGraph,to_workers::Vector{Int64};remote_name = :graph)

    #NOTE: Does not yet support subgraphs.  Aggregate first

    #Create remote channel to store the nodes we want to send
    channel_nodes = RemoteChannel(1)    #we will allocate and send nodes to workers
    channel_indices = RemoteChannel(1)
    channel_master = RemoteChannel(1)   #we will send the master problem to each worker

    n_nodes = getnumnodes(mg)
    n_workers = length(to_workers)
    nodes_per_worker = Int64(floor(n_nodes/n_workers))
    nodes = getnodes(mg)
    node_indices = [getindex(mg,node) for node in nodes]

    link_data = ModelGraphs.get_link_constraint_data(mg)
    n_linkeq_cons = length(link_data.linear_eq_constraints)
    n_linkineq_cons = length(link_data.linear_le_constraints) + length(link_data.linear_ge_constraints) + length(link_data.linear_interval_constraints)
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
    master = getmastermodel(mg)
    println("putting master in channel")
    put!(channel_master, [master])  #put master model in channel

    println("Distributing graph among workers: $to_workers")
    remote_references = []
    #Fill channel with sets of nodes to send
    for (i,worker) in enumerate(to_workers)
        println(worker)
        @spawnat(1, put!(channel_nodes, allocations[i]))
        @spawnat(1, put!(channel_indices, node_indices[i]))
        ref1 = @spawnat worker begin
            Core.eval(Main, Expr(:(=), :master, fetch(channel_master)[1]))
            Core.eval(Main, Expr(:(=), :nodes, take!(channel_nodes)))
            Core.eval(Main, Expr(:(=), :node_indices, take!(channel_indices)))
        end
        wait(ref1)

        ref2 = @spawnat worker Core.eval(Main, Expr(:(=), remote_name, ModelGraphs._create_worker_modelgraph(getfield(Main,:master),getfield(Main,:nodes),getfield(Main,:node_indices),n_nodes,n_linkeq_cons,n_linkineq_cons)))
        push!(remote_references,ref2)
    end

    #NOTE: Linkconstraints keep their indices in new graphs, NOTE: Link constraint row index needs to match on each worker

    return remote_references
end
