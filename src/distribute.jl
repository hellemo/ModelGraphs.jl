
#Distribute a modelgraph among workers.  Each worker should have the same master model.  Each worker will be allocated some of the nodes in the modelgraph
function distribute(mg::ModelGraph,to_workers::Vector{Int64})

    #NOTE: Does not yet support subgraphs.  Aggregate first
    #Create remote channel to store the nodes we want to send
    r_nodes = RemoteChannel(1)
    r_master = RemoteChannel(1)

    n_nodes = getnumnodes(mg)
    n_workers = length(to_workers)
    nodes_per_worker = Int64(floor(n_nodes/n_workers))
    nodes = getnodes(mg)

    allocations = []      #Al4(locate modelnodes onto defined workers
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
        ref = @spawnat(worker, Core.eval(Main, Expr(:(=), :nodes, fetch(r_nodes))))
        # ref2 = @spawnat(worker, Core.eval(Main, Expr(:(=), :master, fetch(r_master))))
        # ref3 = @spawnat(worker,Core.eval(Main, Expr(:(=), :mg, create_modelgraph(getfield(Main,:master),getfield(Main,:nodes)))))
        push!(remote_references,ref)
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
