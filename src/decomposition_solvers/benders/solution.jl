mutable struct BendersSolution <:AbstractSolution
    problemname
    method::Symbol
    solvetime
    objval::Float64
    bestbound::Float64
    gap::Float64
    numiterations::Int64
    termination

    # Iteration Data
    iterval::Array{Float64,1}
    iterbound::Array{Float64,1}
    itertime::Array{Float64,1} # in seconds
    clocktime::Array{Float64,1}
end
