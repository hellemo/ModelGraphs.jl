using ModelGraphs
using Plots
theme(:juno)

graph = ModelGraph()

@linkvariable(graph,z[1:2])
@masterconstraint(graph,z[1] + z[2] <= 2)

#Add nodes to a GraphModel
n1 = add_node!(graph)
n2 = add_node!(graph)

@variable(n1,0 <= x <= 2)
@variable(n1,0 <= y <= 3)
@variable(n1, z >= 0)
@constraint(n1,x+ y + z >= 4)
link_variables!(graph[:z][1],n1[:z])

@variable(n2,x)
@NLnodeconstraint(n2,ref,exp(x) >= 2)
@variable(n2,z >= 0)
@constraint(n2,z + x >= 4)
link_variables!(graph[:z][2],n2[:z])

@linkconstraint(graph,n1[:x] == n2[:x])
@graphobjective(graph,Min,n1[:y] + n2[:x])


p = plotblockmatrix(graph,markersize = 15)
