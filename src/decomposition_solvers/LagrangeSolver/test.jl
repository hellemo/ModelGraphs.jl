using Revise
using JuMP
using GLPK
using AlgebraicGraphs
using LagrangeSolver

m1 = Model(with_optimizer(GLPK.Optimizer))
@variable(m1, xm[i in 1:2],Bin)
@constraint(m1, xm[1] + xm[2] <= 1)
@objective(m1, Max, 16xm[1] + 10xm[2])

m2 = Model(with_optimizer(GLPK.Optimizer))
@variable(m2, xs[i in 1:2], Bin)
@variable(m2, y[i in 1:2], Bin)
@constraint(m2, y[1] + y[2] <= 1)
@constraint(m2, 8xs[1] + 2xs[2] + y[1] + 4y[2] <= 10)
@objective(m2, Max, 4y[2])

m3 = Model(with_optimizer(GLPK.Optimizer))
@variable(m3, x3[i in 1:2], Bin)
@variable(m3, y[i in 1:2], Bin)
@constraint(m3, y[1] + y[2] <= 1)
@constraint(m3, 8x3[1] + 2x3[2] + y[1] + 4y[2] <= 10)
@objective(m3, Max, y[2])

## Model Graph
graph = ModelGraph()

#Add nodes
n1 = add_node!(graph,m1)
n2 = add_node!(graph,m2)
n3 = add_node!(graph,m3)

#link constraints between models
@linkconstraint(graph, [i in 1:2], n1[:xm][i] == n2[:xs][i])
@linkconstraint(graph, n3[:x3][1] + n1[:xm][1] + n2[:xs][1] <= 2)
@linkconstraint(graph, n1[:xm][2] <= n3[:y][1])

lmodel = LagrangeModel(graph)
LagrangeSolver.solve(lmodel)

#print solution

# #verify with full problem
# glpk = with_optimizer(GLPK.Optimizer)
# m,ref_map = create_jump_graph_model(graph)
# optimize!(m,glpk)
