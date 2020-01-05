using JuMP
using ModelGraphs
using Ipopt

function nl_model()
    m = Model()
    @variable(m,x >= rand())
    @variable(m,y >= 1)
    @constraint(m,x + y <= 5)
    @NLconstraint(m,exp(x) >= 2)
    @objective(m,Min,x + y)
    return m
end

#the top level graph
graph = ModelGraph()

#System 1
subgraph1 = ModelGraph()
n1 = add_node!(subgraph1,nl_model())
n2 = add_node!(subgraph1,nl_model())
@linkconstraint(subgraph1,n1[:x] == n2[:x])  #linkconstraint is local to graph1

#System 2
subgraph2 = ModelGraph()
n3 = add_node!(subgraph2,nl_model())
n4 = add_node!(subgraph2,nl_model())
@linkconstraint(subgraph2,n3[:x] == n4[:x])

#Top level links
#check which link constraints I can get here (do I need a getalllinkconstraints function?)
add_subgraph!(graph,subgraph1)
add_subgraph!(graph,subgraph2)
@linkconstraint(graph,n1[:x] == n3[:x])
@linkconstraint(graph,n1[:x] <= n3[:x])

optimize!(graph,with_optimizer(Ipopt.Optimizer))

println("n1[:x]= ",nodevalue(n1[:x]))
println("n2[:x]= ",nodevalue(n2[:x]))
println("n3[:x]= ",nodevalue(n3[:x]))
println("n4[:x]= ",nodevalue(n4[:x]))

println("objective = ", objective_value(graph))
