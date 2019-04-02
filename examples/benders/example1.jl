using JuMP
using GLPKMathProgInterface
using AlgebraicGraphs

#Create the ModelGraph
graph = ModelGraph(with_optimizer = BendersSolver(lp_solver = GLPKSolverLP()))
@variable(graph,y >= 0) #this is a graph variable

#Create the subproblem
subproblem = Model(solver = GLPKSolverLP())
node = addnode!(graph)
setmodel(node,sub_problem)
@variable(subproblem,x[1:2]>=0)
@variable(subproblem,y>=0)

# Indicate a subproblem variable is linked with the graph variable
#2 ways to do this
linkvariable!(y,graph[:y])   #annotate the y variable in the subproblem. It gets promoted to a graph variable.  y would become a parameter in a Benders algorithm
@constraint(subproblem,c1,2*x[1]-x[2]+3*graph[:y] >= 4)

#OR TODO @constraint(subproblem,2*x[1] - x[2] + 3*graph[:y] >= 4)

#Subproblem constraint and objective
@constraint(subproblem,x[1]+2x[2]+y>=3)
@objective(subproblem,Min,2*x[1]+3*x[2])


solve(graph, max_iterations=5)  #Solve with Benders Decomposition

#Query
getgraphvariables(subproblem)  #Returns variables in the sub-problem which are shared with the master level graph

#Change graph objective function
@graphobjective(graph,Min,y + getobjective(node))  #Set the graph objective using the sum of subproblems and the graph variables.
solve(graph)

#Graph Transformations
# BendersSolver works with a SharedVariableGraph.  (No Link Constraints.  Only Shared Variables.)
# LagrangeSolver works with SharedConstraintGraph.  (No shared variables, only shared constraints.)
# PipsSolver works with any ModelGraph.  (Both shared variables and shared constraints.)
