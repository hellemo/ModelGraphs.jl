using Test
#Test basic model functionality
#println("Testing Basic Model Functions")


println("Test adding nonlinear objectives")
@test include("add_NL_objectives.jl")
@test include("nl_problem.jl")

#Special Test Problems
println("Running special test cases")
@test include("test_problems/test_stochastic_pid_problem.jl")
