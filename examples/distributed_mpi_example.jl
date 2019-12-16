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

remote_references = distribute(graph,workers,:graph)  #create the variable graph on each worker


@mpi_do manager begin
    # get local mg
    pipsnlp_solve(graph)
end
