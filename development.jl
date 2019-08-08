using Revise
using Pkg
Pkg.activate(".")

using AlgebraicGraphs

mg = ModelGraph()

n1 = add_node!(mg)
@variable(n1,x[1:5])
@constraint(n1,x[1] + x[2] <= 3)

@variable(mg,z)
