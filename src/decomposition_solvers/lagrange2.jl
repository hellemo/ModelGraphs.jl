mutable struct LagrangeSolverData
    max_iterations::Int
    update_method::Symbol
    epsilon::Float64
    time_limit::Float64
    alpha::Float64
    heuristic::Function  #e.g. fix binaries, fix master problem solution
    delta::Float64
    maxnoimprove::Float64
    cpbound::Float64  #cutting plane upperbound if using cutting plane update method
    node_solver::Union{Nothing,JuMP.OptimizerFactory}  #ClpSolver
    cutting_plane_solver::Union{Nothing,JuMP.OptimizerFactory}

    function LagrangeSolverData()
        data = new()
        data.max_iterations = 10
        data.update_method = SubGradientUpdate()
        data.epsilon = 0.001
        data.timelimit = 3600
        data.alpha = 2
        data.heurstic = FixBinaries()
        data.delta = 0.5
        data.cutting_plane_solver = nothing
        data.node_solver = nothing  #ClpSolver()
        data.solution = LagrangeSolution()
        return data
    end
end

mutable struct LagrangeModel
    subproblems::Vector{JuMP.Model}
    shared_variables::Vector{JuMP.AbstractJuMPScalar}
    link_eq_constraints::Vector{GraphConstraintRef}
    link_ineq_constraints::Vector{GraphConstraintRef}
end

function lagrangesolve(graph::ModelGraph,solver = LagrangeSolver())

    ### INITIALIZATION ###
    lagrange_initialize!(graph,solver.solver_data)

    #switch objective signs
    normalized = _normalize_graph(graph)

    #n = getattribute(graph,:normalized)  #need to track which node indices had their objectives switched
    #initialize multipliers with lp relaxation
    if initialmultipliers == :lp_relaxation
        solve_lp_relaxation!(graph)
    end

    starttime = time()

    #create solution object
    s = Solution(method=:dual_decomposition)

    lambda = graph[:lambda]
    x = graph[:x]

    λ = getattribute(graph,:λ)[end]  #latest multipliers
    x = getattribute(graph,:x)[end]  #latest primal solution

    residual = getattribute(graph,:res)[end]
    n_multipliers = getattribute(graph,:numlinks)
    #nodes = [node for node in getnodes(graph)]
    setattribute(graph,:α, [1.0*α])
    iterval = 0

    ### ITERATIONS ###
    for iter in 1:max_iterations
        variant = iter == 1 ? :default : update_method # Use default method in the first iteration

        iterstart = time()
        # Solve subproblems
        Zk = 0  #objective value
        #NOTE: We can parallelize here
        for node in getnodes(graph)
           (x,Zkn) = solvenode(node,λ,x,variant) #x gets updated by solvenode.
           Zk += Zkn        #add objective values
        end
        setattribute(graph, :steptaken, true)

        # If no improvement, increase no improvement counter
        if Zk < getattribute(graph, :Zk)[end]
            #graph[:noimprove] = graph[:noimprove] + 1
            setattribute(graph, :noimprove, getattribute(graph, :noimprove) + 1)
        end
        # If too many iterations without improvement, decrease :α
        if getattribute(graph, :noimprove) >= getattribute(graph, :maxnoimprove)
            setattribute(graph, :noimprove, 0)
            getattribute(graph, :α)[end] *= getattribute(graph, :δ)
        end
        # Save info
        push!(getattribute(graph, :Zk),Zk)
        push!(getattribute(graph, :x),x)
        α = getattribute(graph, :α)[end]
        push!(getattribute(graph, :α), α)

        # Update residuals
        res = x[:,1] - x[:,2]
        push!(getattribute(graph, :res),res)

        # Save iteration data
        itertime = time() - iterstart
        tstamp = time() - starttime
        saveiteration(s,tstamp,[n*iterval,n*Zk,itertime,tstamp],n)
        printiterationsummary(s,singleline=false)

        # Check convergence
        if norm(res) < ϵ
            s.termination = "Optimal"
            return s
        end

        # Check time limit
        if tstamp > timelimit
            s.termination = "Time Limit"
            return s
        end

        # Update multipliers
        println("α = $α")
        (λ, iterval) = updatemultipliers(graph,λ,res,update_method,lagrangeheuristic)
        push!(getattribute(graph, :λ), λ)
        # Update iteration time
        s.itertime[end] = time() - iterstart
    end
    s.termination = "Max Iterations"
    return s
end
