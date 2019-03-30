using JuMP
using GLPKMathProgInterface
using Plasmo

##Place MP and SP into ModelGraph
# mp = Model(solver = GLPKSolverLP())
# @variable(mp,y>=0)
# @objective(mp,Min,2y)

graph = ModelGraph(with_optimizer = BendersSolver(lp_solver = GLPKSolverLP()))
@graphvariable(graph,y >= 0)

subproblem = Model(solver = GLPKSolverLP())
@variable(subproblem,x[1:2]>=0)
@variable(subproblem,y>=0)
@constraint(subproblem,c1,2x[1]-x[2]+3y>=4)
@constraint(subproblem,x[1]+2x[2]+y>=3)
@objective(subproblem,Min,2*x[1]+3x[2])

node = addnode!(graph)
setmodel(node,sub_problem)

set_shared_variable(graph[:y],node[:y])  #annotate the y variable in the subproblem. It gets promoted to a graph variable.  y would become a parameter in a Benders algorithm

solve(graph, max_iterations=5)

#We can also do things like:
@constraint(subproblem,x[1]+2*x[2]+graph[:y] >= 3)  #This would create a graph constraint

getgraphvariables(subproblem)  #Returns variables in the sub-problem which are shared with the master level graph

#Graph Transformations
# BendersSolver works with a SharedVariableGraph.  (No Link Constraints.  Only Shared Variables.)
# LagrangeSolver works with SharedConstraintGraph.  (No shared variables, only shared constraints.)
# PipsSolver works with any ModelGraph.  (Both shared variables and shared constraints.)
