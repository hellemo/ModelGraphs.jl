#Distribute a modelgraph among workers.  Each worker should have the same master model.  Each worker will be allocated some of the nodes in the modelgraph
function distribute(mg::ModelGraph,to_workers::Vector{Int64};remote_name = :graph)

    #NOTE: Does not yet support subgraphs.  Aggregate first
    #Create remote channel to store the nodes we want to send
    r_nodes = RemoteChannel(1)    #we will allocate and send nodes to workers
    r_master = RemoteChannel(1)   #we will send the master problem to each worker

    n_nodes = getnumnodes(mg)
    n_workers = length(to_workers)
    nodes_per_worker = Int64(floor(n_nodes/n_workers))
    nodes = getnodes(mg)

    link_data = get_link_constraint_data(mg)

    n_linkeq_cons = length(link_data.linear_eq_constraints)
    n_linkineq_cons = length(link_data.linear_le_constraints) + length(link_data.linear_ge_constraints) + length(link_data.linear_interval_constraints)
    n_link_cons = n_linkeq_cons + n_linkineq_cons

    master = getmastermodel(mg)
    put!(r_master, master)  #put master model in channel

    #Allocate modelnodes onto provided workers
    allocations = []
    j = 1
    while  j < n_nodes
        if j + nodes_per_worker > n_nodes
            push!(allocations,nodes[j:end])
        else
            push!(allocations,nodes[j:j+nodes_per_worker - 1])
        end
        j += nodes_per_worker
    end

    println("Distributing graph among workers: $to_workers")
    remote_references = []
    #Fill channel with sets of nodes to send
    for (i,worker) in enumerate(to_workers)
        @spawnat(1, put!(r_nodes, allocations[i]))
        #@spawnat worker begin
        #Core.eval(Main, Expr(:(=), :nodes, fetch(r_nodes)))
        #Core.eval(Main, Expr(:(=), :master, fetch(r_master))
        #end
        #ref = @spawnat worker Core.eval(Main, Expr(:(=), $remote_name, ModelGraphs._create_worker_modelgraph(getfield(Main,:master),getfield(Main,:nodes),n_linkeq_cons,n_linkineq_cons)))
        @spawnat(worker, Core.eval(Main, Expr(:(=), :nodes, fetch(r_nodes))))
        @spawnat(worker, Core.eval(Main, Expr(:(=), :master, fetch(r_master))))

        #WAIT here
        ref3 = @spawnat(worker, Core.eval(Main, Expr(:(=), remote_name, ModelGraphs._create_worker_modelgraph(getfield(Main,:master),getfield(Main,:nodes),n_linkeq_cons,n_linkineq_cons))))
        push!(remote_references,ref3)
    end

    #NOTE: Linkconstraints keep their indices in new graphs, NOTE: Link constraint row index needs to match on each worker

    return remote_references
end



# function send_graph_data(manager::MPI.MPIManager,graph::ModelGraph)
#
#     julia_workers = collect(values(manager.mpi2j))
#     r = RemoteChannel(1)
#
#     @spawnat(1, put!(r, [graph]))
#
#     @sync for to in julia_workers
#         @spawnat(to, Core.eval(Main, Expr(:(=), :graph, fetch(r)[1])))
#     end
# end