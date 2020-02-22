using ModelGraphs
using Ipopt

graph = ModelGraph()
optimizer = Ipopt.Optimizer

#m1 = JuMP.Model()
n1 = add_node!(graph)
@variable(n1,0 <= x <= 2)
@variable(n1,0 <= y <= 3)
@constraint(n1,x+y <= 4)
@objective(n1,Min,2*x)

#m2 = Model()
n2 = add_node!(graph)
@variable(n2,x)
@NLnodeconstraint(n2,exp(x) >= 2)

#Link constraints take the same expressions as the JuMP @constraint macro
@linkconstraint(graph,n1[:x] == n2[:x])

agg_model,reference_map = aggregate(graph)
optimize!(agg_model,optimizer)

#Use the reference map to look up values
println("n1[:x] = ",value(reference_map[n1[:x]]))
println("n2[:x] = ",value(reference_map[n2[:x]]))
