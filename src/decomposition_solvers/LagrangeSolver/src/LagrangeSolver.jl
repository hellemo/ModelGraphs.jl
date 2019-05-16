module LagrangeSolver

using AlgebraicGraphs
using JuMP
using MathOptInterface
const MOI = MathOptInterface
using SparseArrays
using LinearAlgebra

export LagrangeModel, dual_decomposition_solve, solve

include("utils.jl")

include("solution.jl")

include("lagrange_model.jl")

include("multiplier_updates.jl")

end
