using JuMP
using ModelGraphs
using Ipopt

graph = ModelGraph()

#Add nodes to a GraphModel
n1 = add_node!(graph)
n2 = add_node!(graph)
n3 = add_node!(graph)
n4 = add_node!(graph)
#Add edges between the nodes

#Set a model on node 1
m1 = Model()
@variable(m1,0 <= x <= 2)
@variable(m1,0 <= y <= 3)
@NLconstraint(m1,x^3+y <= 4)
@objective(m1,Min,x)

#Set a model on node 2
m2 = Model()
vals = collect(1:5)
grid = 1:3
@variable(m2,x >= 1)
@variable(m2,0 <= y <= 5)
@variable(m2,z[1:5] >= 0)
@variable(m2,a[vals,grid] >=0 )
@NLconstraint(m2,exp(x)+y <= 7)
@objective(m2,Min,x)

m3 = Model()
@variable(m3,x[1:5])

m4 = Model()
@variable(m4,x <= 1)


#Set models on nodes and edges
set_model(n1,m1)     #set m1 to node 1.  Updates reference on m1
set_model(n2,m2)
set_model(n3,m3)
set_model(n4,m4)

#Link constraints take the same expressions as the JuMP @constraint macro
@linkconstraint(graph,n4[:x] == n1[:x])
@linkconstraint(graph,[t = 1:5],n4[:x] == n2[:z][t])
@linkconstraint(graph,[i = 1:5],n3[:x][i] == n1[:x])
@linkconstraint(graph,[j = 1:5,i = 1:3],n2[:a][j,i] == n4[:x])
@linkconstraint(graph,[i = 1:3],n1[:x] + n2[:z][i] + n3[:x][i] + n4[:x] >= 0)

ipopt = Ipopt.Optimizer
optimize!(graph,ipopt)

#Query solution
@test nodevalue(n2[:x]) ≈ 1
