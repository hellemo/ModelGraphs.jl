using ModelGraphs
using Ipopt

graph = ModelGraph()
optimizer = with_optimizer(Ipopt.Optimizer)

#Master node
@variable(graph,z[1:2])
@constraint(graph,z[1] + z[2] <= 2) #constraint goes to a master node

#Add nodes to a ModelGraph
n1 = add_node!(graph)
n2 = add_node!(graph)

#Node 1 Model
@variable(n1,0 <= x <= 2)
@variable(n1,0 <= y <= 3)
@variable(n1, z >= 0)
@constraint(n1,x+y+z >= 4)
@link graph graph[:z][1] n1[:z]

#Node 2 Model
@variable(n2,x)
@NLnodeconstraint(n2,ref,exp(x) >= 2)
@variable(n2,z >= 0)
@constraint(n2,z + x >= 4)
@link graph graph[:z][2] n2[:z]

#Link constraints take the same expressions as the JuMP @constraint macro
@linkconstraint(graph,n1[:x] == n2[:x])

#Graph objective
@graphobjective(graph,Min,n1[:y] + n2[:x])

#Optimize with Ipopt.  (This will aggregate the graph into an Ipopt Model)
optimize!(graph,optimizer)

println("n1[:z]= ",nodevalue(n1[:z]))
println("n2[:z]= ",nodevalue(n2[:z]))
println("n1[:x]= ",nodevalue(n1[:x]))
println("n1[:y]= ",nodevalue(n1[:y]))
println("n2[:x]= ",nodevalue(n2[:x]))
println("objective = ", objective_value(graph))
