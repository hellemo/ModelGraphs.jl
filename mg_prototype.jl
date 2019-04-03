using Pkg
Pkg.activate("./")

using Revise
using JuMP
using AlgebraicGraphs

graph = ModelGraph()
@variable(graph,z)         # This is a graph variable
@constraint(graph,z <= 4)  # This is a graph constraint

#Node 1
m1 = Model()
@variable(m1,x)
n1 = addnode!(graph,m1)

#Node 2
m2 = Model()
@variable(m2,y)
n2 = addnode!(graph,m2)

@linkconstraint(graph,n1[:x] + n1[:y] <= 4)  #This a link constraint
