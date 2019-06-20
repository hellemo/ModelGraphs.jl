using JuMP
using AlgebraicGraphs
using Ipopt

const AG = AlgebraicGraphs

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
@constraint(m1,ref,x+y <= 4)
@objective(m1,Min,x)

#Set a model on node 2
m2 = Model()
vals = collect(1:5)
grid = 1:3
@variable(m2,x >= 1)
@variable(m2,0 <= y <= 5)
@variable(m2,z[1:5] >= 0)
@variable(m2,a[vals,grid] >=0 )
@NLconstraint(m2,nl_con_ref,exp(x)+y <= 7)
@objective(m2,Min,x)

m3 = Model()
@variable(m3,x[1:5])

m4 = Model()
@variable(m4,x <= 1)


#Set models on nodes and edges
setmodel(n1,m1)     #set m1 to node 1.  Updates reference on m1
setmodel(n2,m2)
setmodel(n3,m3)
setmodel(n4,m4)

#Link constraints take the same expressions as the JuMP @constraint macro
@linkconstraint(graph,n4[:x] == n1[:x])
@linkconstraint(graph,[t = 1:5],n4[:x] == n2[:z][t])
@linkconstraint(graph,[i = 1:5],n3[:x][i] == n1[:x])
@linkconstraint(graph,[j = 1:5,i = 1:3],n2[:a][j,i] == n4[:x])
@linkconstraint(graph,[i = 1:3],n1[:x] + n2[:z][i] + n3[:x][i] + n4[:x] >= 0)

# Now Aggregate 4 models together
nodes = collect(getnodes(graph))
edges = collect(getedges(graph))
aggregated_model,reference_map = AG.create_aggregate_model(graph,nodes,edges)
optimize!(aggregated_model,with_optimizer(Ipopt.Optimizer))

#Use reference map to look up aggregated variable values
for var in JuMP.all_variables(m1)
    println(value(reference_map[var]))
end

#Do the same thing to look up m2 values in aggregated model
for var in JuMP.all_variables(m2)
    println(var," = ",value(reference_map[var]))
end

constraint_types = JuMP.list_of_constraint_types(m1)
for (func,set) in constraint_types
    constraint_refs = JuMP.all_constraints(m1,func,set)
    for ref in constraint_refs
        aggregate_model_ref = reference_map[ref]
        constraint = constraint_object(aggregate_model_ref)
        println(constraint)
    end
end

#an aggregated model has an associated graph
aggregate_graph = aggregated_model.ext[:Graph]

n1 = getnode(aggregate_graph,1)
n2 = getnode(aggregate_graph,2)
n1.obj_dict
n2.obj_dict
