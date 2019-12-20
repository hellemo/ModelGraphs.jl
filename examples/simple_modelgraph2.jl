using ModelGraphs
#using Ipopt

graph = ModelGraph()
#optimizer = with_optimizer(Ipopt.Optimizer)

# @linkvariable(graph,z[1:2])
# @masterconstraint(graph,z[1] + z[2] <= 2)

#Add nodes to a GraphModel
n1 = add_node!(graph)
n2 = add_node!(graph)
n3 = add_node!(graph)
n4 = add_node!(graph)

@variable(n1,0 <= x <= 2)
@variable(n1,0 <= y <= 3)
@variable(n1, z >= 0)
@constraint(n1,x+y+z >= 4)
@objective(n1,Min,y)

@variable(n2,x >= 0)
@NLnodeconstraint(n2,ref,exp(x) >= 2)
@variable(n2,z >= 0)
@constraint(n2,z + x >= 4)
@objective(n2,Min,x)

@variable(n3,x[1:5] >= 0)
@NLnodeconstraint(n3,ref,exp(x[3]) >= 5)
@constraint(n3,sum(x[i] for i = 1:5) == 10)
@objective(n3,Min,x[1] + x[2] + x[3])

@variable(n4,x[1:5])
@constraint(n4,sum(x[i] for i = 1:5) >= 10)
@NLnodeconstraint(n4,ref,exp(x[2]) >= 4)
@objective(n4,Min,x[2])

#Link constraints take the same expressions as the JuMP @constraint macro
@linkconstraint(graph,n1[:x] == n2[:x])
@linkconstraint(graph,n2[:x] == n3[:x][3])
@linkconstraint(graph,n3[:x][1] <= n4[:x][1])
@linkconstraint(graph,n1[:x] + n2[:x] <= n4[:x][5])
