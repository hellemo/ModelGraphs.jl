using ModelGraphs
using JuMP
using Ipopt
using MathOptInterface
const MOI = MathOptInterface

println("Testing addition of nonlinear objective functions")

graph = ModelGraph()

m1 = Model()
@variable(m1,x[1:2] <= 2)
set_start_value(x[1],2)
set_start_value(x[2],1)
@NLobjective(m1,Max,x[1]^2 + x[2]^2)
n1 = add_node!(graph,m1)

m2 = Model()
@variable(m2,x[1:2] >= 0 )
setvalue(x[1],2)
setvalue(x[2],2)
@NLobjective(m2,Min,x[1]^3 + x[2]^2)
n2 = add_node!(graph,m2)

aggmodel,ref_map = aggregate(graph)
optimize!(aggmodel,with_optimizer(Ipopt.Optimizer))

d = JuMP.NLPEvaluator(aggmodel)
MOI.initialize(d,[:ExprGraph])

true
