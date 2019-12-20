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

#Get references to the graph on each worker
remote_references = distribute(mg, workers(), remote_name = :graph)
r1 = remote_references[1]
r2 = remote_references[2]
