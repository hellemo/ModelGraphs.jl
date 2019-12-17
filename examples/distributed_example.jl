# to import MPIManager
#using MPIClusterManagers

# need to also import Distributed to use addprocs()
using Distributed

if nprocs() == 1
    addprocs(2)
end
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using ModelGraphs

include("simple_modelgraph.jl")


# for worker in workers()
#     ref1 = @spawnat(worker, Core.eval(Main, Expr(:(=), :master, 2)))
#     wait(ref1)
# end

remote_references = distribute(graph, workers(), remote_name = :graph)

r1 = remote_references[1]
r2 = remote_references[2]


# # #Test that create worker graph is working
# mg = graph
# master = getmastermodel(mg)
# nodes = getnodes(mg)
# link_data = ModelGraphs.get_link_constraint_data(mg)
# n_linkeq_cons = length(link_data.linear_eq_constraints)
# n_linkineq_cons = length(link_data.linear_le_constraints) + length(link_data.linear_ge_constraints) + length(link_data.linear_interval_constraints)
#
# f1 = @spawnat 2 ModelGraphs._create_worker_modelgraph(getfield(Main,:master),getfield(Main,:nodes),getfield(Main,:node_indices),2,n_linkeq_cons,n_linkineq_cons)
# f2 = @spawnat 3 ModelGraphs._create_worker_modelgraph(getfield(Main,:master),getfield(Main,:nodes),getfield(Main,:node_indices),2,n_linkeq_cons,n_linkineq_cons)
#
# mg_new = ModelGraphs._create_worker_modelgraph(master,nodes,[1,2],2,n_linkeq_cons,n_linkineq_cons)


# specify, number of mpi workers, launch cmd, etc.
# manager=MPIManager(np=2)
# # start mpi workers and add them as julia workers too.
# addprocs(manager)
# @everywhere using ModelGraphs
#
#
# #Create ModelGraph
# mg = ModelGraph()
#
# #Create model
# #
# #
# #
#
#
# refs = distribute!(mg,procs(manager))   #Create a remote reference to each new graph
#
# @mpi_do manager begin
#     # get local mg
#
#
#     pipsnlpsolve(mg)
# end
