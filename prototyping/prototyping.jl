using Pkg
Pkg.activate("./")

using Revise
using JuMP
using AlgebraicGraphs

lm = AlgebraicGraphs.LinkModel()

@variable(lm,z)         # This is a graph variable
@constraint(lm,z <= 4)  # This is a graph constraint

m1 = Model()
@variable(m1,x)

m2 = Model()
@variable(m2,y)

#This a link constraint
@constraint(lm,x + y <= 4)
