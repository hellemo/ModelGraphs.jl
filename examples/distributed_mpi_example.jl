# to import MPIManager
using MPIClusterManagers

# need to also import Distributed to use addprocs()
using Distributed

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using ModelGraphs
#@everywhere using ModelGraphMPISolvers

include("simple_modelgraph1.jl")

# specify, number of mpi workers, launch cmd, etc.
if !(isdefined(Main,:manager))
    manager=MPIManager(np=2)
    # start mpi workers and add them as julia workers too.
    addprocs(manager)
end
#NOTE: This seems to have no value the first time I run this.
julia_workers = collect(values(manager.mpi2j))

#Distribute the graph to workers
remote_references = distribute(graph,julia_workers,remote_name = :graph)  #create the variable graph on each worker

@mpi_do manager begin
    using MPI
    println(graph)
end

rank_zero = manager.mpi2j[0] #julia process representing rank 0
