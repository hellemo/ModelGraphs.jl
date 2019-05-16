module Lagrange

using AlgebraicGraphs
using JuMP
using MathOptInterface
const MOI = MathOptInterface
using SparseArrays

export LagrangeModel, lagrange_solve

include("utils.jl")

include("solution.jl")

include("dual_decomposition.jl")

include("multiplier_updates.jl")
