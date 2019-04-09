using Pkg
Pkg.activate("./")

using Revise
using JuMP
using AlgebraicGraphs
using Ipopt

graph = ModelGraph()
#Node 1
m1 = Model()
@variable(m1,x)
n1 = add_node!(graph,m1)

#Node 2
m2 = Model()
@variable(m2,y)
n2 = add_node!(graph,m2)

@linkconstraint(graph,n1[:x] + n2[:y] <= 4)  #This a link constraint

#optimizer = with_optimizer(Ipopt.Optimizer)
