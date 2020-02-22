using Test
#Test basic model functionality
#println("Testing Basic Model Functions")


@testset "Test adding nonlinear objectives" begin
    include("add_NL_objectives.jl")
    include("nl_problem.jl")
end

@testset "Running special test cases" begin
    @test_broken include("test_problems/test_stochastic_pid_problem.jl")
end
