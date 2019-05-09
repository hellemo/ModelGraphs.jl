#NOTE: This file contains proposed (not yet implemented syntax)

using JuMP
using Ipopt
using Plasmo
using PlasmoSolverInterface #gets PipsSolver

function get_electricity_model(demand)
    m = Model()
    #amount of electricity produced
    @variable(m, 0<=prod<=10, start=5)
    #amount of electricity purchased or sold
    @variable(m, input, start = 2)
    #amount of gas purchased
    @variable(m, gas_purchased, start = 2)
    @constraint(m, gas_purchased >= prod)
    @constraint(m, prod + input == demand)
    return m
end

Ns = 15
demand = rand(Ns)*10
graph = ModelGraph()
@graphvariable(graph,0<=gas_purchased<=8, start = 2)
@objective(master,Min,gas_purchased)

for j in 1:Ns
    scenm = get_electricity_model(demand[j])
    node = addnode!(graph,scenm)
    linkvariable!(graph[:gas_purchased],node[:gas_purchased])
    #connect children and parent variables
    #@linkconstraint(graph, master[:gas_purchased] == scenm[j][:gas_purchased])
    #Create child objective
    @objective(scenm[j],Min,1/Ns*(scenm[j][:prod] + 3*scenm[j][:input]))

end

#create a link constraint between the subproblems (PIPS-NLP supports this kind of constraint)
@linkconstraint(graph, (1/Ns)*sum(scenm[s][:prod] for s in 1:Ns) == 8)
@graphobjective(graph,Min,graph[:gas_purchased] + sum(node_objectives(graph)))

solver = PipsSolver(n_workers = 2)
solve(graph,solver)
