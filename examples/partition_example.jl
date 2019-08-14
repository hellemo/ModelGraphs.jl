using JuMP
using AlgebraicGraphs  #ModelGraphs?
using Ipopt
using KaHyPar
using SparseArrays

graph = ModelGraph()

#setsolver(graph,Ipopt.IpoptSolver())

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
@constraint(m1,x+y <= 4)
@objective(m1,Min,x)

#Set a model on node 2
m2 = Model()
@variable(m2,x >= 1)
@variable(m2,0 <= y <= 5)
@NLconstraint(m2,exp(x)+y <= 7)
@objective(m2,Min,x)

m3 = Model()
@variable(m3,x >= 0)

m4 = Model()
@variable(m4,0 <= x <= 1)


#Set models on nodes and edges
set_model(n1,m1)     #set m1 to node 1.  Updates reference on m1
set_model(n2,m2)
set_model(n3,m3)
set_model(n4,m4)

ipopt = with_optimizer(Ipopt.Optimizer)

#Link constraints take the same expressions as the JuMP @constraint macro
@linkconstraint(graph,n4[:x] == n1[:x])
@linkconstraint(graph,n1[:y] + n2[:y] + n3[:x] <= 2 )

hypergraph = gethypergraph(graph)
A = sparse(hypergraph)
partition1 = KaHyPar.partition(A,2,configuration = :edge_cut)
partition2 = KaHyPar.partition(A,2,configuration = :connectivity)

optimize!(graph,ipopt)

println("n1[:x]= ",nodevalue(n1[:x]))
println("n1[:y]= ",nodevalue(n1[:y]))
println("n4[:x]= ",nodevalue(n4[:x]))
