using JuMP
using GLPKMathProgInterface
using Plasmo

##Place MP and SP into ModelGraph
# mp = Model(solver = GLPKSolverLP())
# @variable(mp,y>=0)
# @objective(mp,Min,2y)

graph = ModelGraph(with_optimizer = BendersSolver(lp_solver = GLPKSolverLP()))
@variable(graph,y >= 0) #graph variable

subproblem = Model(solver = GLPKSolverLP())

node = addnode!(graph)
setmodel(node,sub_problem)

@variable(subproblem,x[1:2]>=0)
@variable(subproblem,y>=0)
linkvariable!(y,graph[:y])
@constraint(subproblem,c1,2x[1]-x[2]+3y>=4)
#OR
#@constraint(subproblem,2*x[1] - x[2] + 3*graph[:y] >= 4)

@constraint(subproblem,x[1]+2x[2]+y>=3)
@objective(subproblem,Min,2*x[1]+3x[2])



#set_shared_variable(graph[:y],node[:y])  #annotate the y variable in the subproblem. It gets promoted to a graph variable.  y would become a parameter in a Benders algorithm

solve(graph, max_iterations=5)

#We can also do things like:
@constraint(subproblem,x[1]+2*x[2]+graph[:y] >= 3)  #This would create a graph constraint

getgraphvariables(subproblem)  #Returns variables in the sub-problem which are shared with the master level graph

@graphobjective(grap,Min,y + getobjective(node))  #Set the graph objective using the sum of subproblems and the graph variables.

#Graph Transformations
# BendersSolver works with a SharedVariableGraph.  (No Link Constraints.  Only Shared Variables.)
# LagrangeSolver works with SharedConstraintGraph.  (No shared variables, only shared constraints.)
# PipsSolver works with any ModelGraph.  (Both shared variables and shared constraints.)
