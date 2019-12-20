using Distributed

if nprocs() == 1
    addprocs(2)
end

#Setup the worker environments
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using ModelGraphs


include("simple_modelgraph2.jl")

#Distribute the modelgraph to the Julia workers
remote_graphs = distribute(graph,workers(),remote_name = :graph)  #create the variable pipsgraph on each worker

r1 = remote_graphs[1] #reference to the model graph on worker 1
r2 = remote_graphs[2] #reference to the model graph on worker 2

f = @spawnat 2 println(fetch(r1).obj_dict)
