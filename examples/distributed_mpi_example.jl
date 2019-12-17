# to import MPIManager
using MPIClusterManagers

# need to also import Distributed to use addprocs()
using Distributed

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using ModelGraphs

include("simple_modelgraph.jl")

# specify, number of mpi workers, launch cmd, etc.
manager=MPIManager(np=2)
# start mpi workers and add them as julia workers too.
addprocs(manager)

workers = manager.mpi2j
remote_references = distribute(graph,workers,remote_name = :graph)  #create the variable graph on each worker
@mpi_do manager begin
    pipsnlp_solve(graph)
end

rank_zero = manager.mpi2j[0] #julia process representing rank 0
solution = fetch(@spawnat(rank_zero, getfield(Main, :graph)))

# @mpi_do manager begin
#     using MPI
#     comm=MPI.COMM_WORLD
#     println("Hello world, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")
# end
