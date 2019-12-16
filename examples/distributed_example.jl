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

remote_references = distribute(graph,workers())


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
