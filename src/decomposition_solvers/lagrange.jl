mutable struct LagrangeSolverData
    max_iterations::Int
    update_method::Symbol
    epsilon::Float64
    timelimit::Float64
    alpha::Float64
end


"""
  lagrangesolve(graph)

  Solve a ModelGraph object using lagrange decomposition
"""
function lagrangesolve(graph::ModelGraph,solver_data = LagrangeSolverData(),
    max_iterations=10,
    update_method=:subgradient,         # probingsubgradient
    ϵ=0.001,                            # ϵ-convergence tolerance
    timelimit=3600,
    α=2,                                # default subgradient step
    lagrangeheuristic=fixbinaries,      # function to calculate the upper bound
    initialmultipliers=:zero,           # relaxation for LP relaxation
    δ = 0.5,                            # Factor to shrink step when subgradient stuck
    maxnoimprove = 3,
    cpbound=1e6,  # Amount of iterations that no improvement is allowed before shrinking step
    cutting_plane_solver = ClpSolver(),
    node_solver = ClpSolver())

    ### INITIALIZATION ###
    lagrange_initialize!(graph,cutting_plane_solver,node_solver,δ,maxnoimprove,cpbound)

    #switch objective signs
    n = getattribute(graph,:normalized)

    #initialize multipliers with lp relaxation
    if initialmultipliers == :relaxation
        solve_lp_relaxation!(graph)
    end

    starttime = time()

    #create solution object
    s = Solution(method=:dual_decomposition)

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
        res = x[:,1] - x[:,2]   #NOTE: This assumes very simple equality constraints
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

# Initialization
"""
  lg_initialize!(graph::ModelGraph)

  Prepares the graph to apply a lagrange decomposition algorithm
"""
function lg_initialize!(graph::ModelGraph,cutting_plane_solver::JuMP.OptimizerFactory,node_solver::JuMP.OptimizerFactory,
    δ=0.5,
    maxnoimprove=3,
    cpbound=nothing)

    if hasattribute(graph,:preprocessed)
        return true
    end

    n = normalizegraph(graph)

    links = getlinkconstraints(graph)
    nmult = length(links)                   # Number of multipliers

    setattribute(graph, :numlinks, nmult)
    setattribute(graph, :λ, [zeros(nmult)])   # Array{Float64}(nmult)
    setattribute(graph, :x, [zeros(nmult,2)]) # Linking variables values
    setattribute(graph, :res, [zeros(nmult)]) # Residuals
    setattribute(graph, :Zk, [0.0]) # Bounds
    setattribute(graph, :cuts, [])
    setattribute(graph, :δ, δ)
    setattribute(graph, :noimprove, 0)
    setattribute(graph, :maxnoimprove, maxnoimprove)
    setattribute(graph, :explore, [])
    setattribute(graph, :steptaken, false)

    # Create Lagrange Master.  Used for cutting plane multiplier updates
    ms = Model(solver = cutting_plane_solver)
    @variable(ms, η, upperbound=cpbound)  #single cut master
    @variable(ms, λ[1:nmult])
    @objective(ms, Max, η)

    setattribute(graph, :lgmaster, ms)

    # Each node saves its initial objective and sets a solver if it doesn't have one
    for node in getnodes(graph)
        m = getmodel(node)
        # if mn.solver == JuMP.UnsetSolver()
        #     JuMP.setsolver(mn, node_solver)
        # end
        m.ext[:preobj] = JuMP.objective_function(m)
        m.ext[:multmap] = Dict()  #multiplier map
        m.ext[:varmap] = Dict()   #variable map
    end

    # Maps
    # Multiplier map to know which component of λ to take
    # Varmap knows what values to post where

    #IDEA: Take a more graph-based approach here
    for (i,lc) in enumerate(links)

        for j in 1:length(lc.terms.vars)
            var = lc.terms.vars[j]
            #NOTE: This is getting over-written for every variable
            var.m.ext[:multmap][i] = (lc.terms.coeffs[j],lc.terms.vars[j])  ##NOTE: This won't work for .  link[i] --> (coeff[j],var[j]) where j is the last variable in the link constraint
            var.m.ext[:varmap][var] = (i,j)
        end

    end
    setattribute(graph, :preprocessed, true)
end

# Solve a single subproblem (node)
# NOTE: Need to get multiplier, need to store primal result
function solvenode(node::ModelNode,λ,x,variant=:default)
    m = getmodel(node)

    # Restore objective function
    # m.obj = m.ext[:preobj]
    JuMP.set_objective_function(m,m.ext[:preobj])
    m.ext[:lgobj] = m.ext[:preobj]

    # Add dualized part to objective function

    for k in keys(m.ext[:multmap])
        coef = m.ext[:multmap][k][1]
        var = m.ext[:multmap][k][2]
        m.ext[:lgobj] += λ[k]*coef*var #Add to objective

        obj = JuMP.objective_function(m)
        obj += λ[k]*coef*var

        #Add ADMM piece to objective function
        if variant == :ADMM
            j = 3 - m.ext[:varmap][var][2]
            obj += 1/2*(coef*var - coef*x[k,j])^2
        end
        JuMP.set_objective_function(m,obj)

    end
    # Solve sub-problem
    optimize!(m)

    # Pass output
    # This updates the solution on other nodes
    for variable in keys(m.ext[:varmap])
        val = JuMP.value(variable)
        x[m.ext[:varmap][variable]...] = val  #store variable values in x.  varmap[variable]... is the (i,j) pair of coordinates for link
    end
    objval = JuMP.value(m.ext[:lgobj])
    #store node objective value
    setattribute(node, :objective, objval)

    #setattribute(node, :solvetime, getsolvetime(m))

    return x, objval
end

#NOTE: functions I might write here
function dual_decomposition_objective
end

function admm_objective
end

# Multiplier Initialization
function solve_lp_relaxation!(graph::ModelGraph,optimizer::JuMP.OptimizerFactory)
    #if graph.jump_model == nothing
    if !hasattribute(graph,:m_relaxation)
        #setattribute(graph, :mflat, create_jump_graph_model(graph))
        relaxed_model = create_jump_graph_model(graph,with_optimizer = optimizer)
        #Relax integrality and binary constraints
        for variable in JuMP.all_variables(relaxed_model)
            if JuMP.is_binary(variable)
                JuMP.unset_binary(variable)
            end
            if JuMP.is_integer(variable)
                JuMP.unset_integer(variable)
            end
        end
        setattribute(graph,relaxed_model)
    end

    n = getattribute(graph , :normalized)
    n_multipliers = getattribute(graph , :numlinks)
    relaxed_model = getattribute(graph , :m_relaxation)

    JuMP.optimize!(relaxed_model)

    #initialize multipliers
    #NOTE: Shouldn't n depend on the node that flipped objective signs?
    #get duals for the linkconstraints
    link_duals = [JuMP.dual(link) for link in relaxed_model.ext[:Graph].linkconstraints]  #grab the link constraints from the aggregated model
    #duals = n*mf.linconstrDuals[end-nmult+1:end]
    #getattribute(graph , :λ)[end] = n*mf.linconstrDuals[end-nmult+1:end]

    #get the last element because it's a multiplier history?
    setattribute(graph , :λ)[end] = link_duals

    return JuMP.objective_value(relaxed_model)
end


function updatemultipliers(graph::ModelGraph,λ::Vector,residual::Vector,method::Symbol,lagrangeheuristic=nothing)
    if method == :subgradient
        subgradient(graph,λ,residual,lagrangeheuristic)
    elseif method == :probingsubgradient
        probingsubgradient(graph,λ,residual,lagrangeheuristic)
    elseif method == :cuttingplanes
        cuttingplanes(graph,λ,residual)
    elseif  method == :interactive
        interactive(graph,λ,residual,lagrangeheuristic)
    end
end


# Update functions
function subgradient(graph::ModelGraph,λ,residual,lagrangeheuristic)
    α = getattribute(graph , :α)[end]
    n = getattribute(graph , :normalized)

    #bound
    bound = n*lagrangeheuristic(graph)


    Zk = getattribute(graph , :Zk)[end]

    step = α*abs(Zk-bound)/(norm(residual)^2)

    λ += step*res
    return λ,bound
end

function αeval(αv,graph,bound)
    xv = deepcopy(getattribute(graph , :x)[end])
    res = getattribute(graph , :res)[end]
    Zk = getattribute(graph , :Zk)[end]
    λ = getattribute(graph , :λ)[end]
    nodes = [node for node in getnodes(graph)]
    step = abs(Zk-bound)/(norm(res)^2)
    zk = 0
    for node in nodes
        (xv,Zkn) = solvenode(node,λ+αv*step*res,xv,:default)
        zk += Zkn
    end
    return zk
end

function αexplore(graph,bound)
    df = getattribute(graph , :explore)
    n = getattribute(graph , :normalized)
    z = Float64[]
    for α in 0:0.1:2
        push!(z,n*αeval(α,graph,bound))
    end
    push!(df,z)
end

function probingsubgradient(graph,λ,res,lagrangeheuristic,α=getattribute(graph , :α)[end],Δ=0.01;exhaustive=false)
    res = getattribute(graph , :res)[end]
    Zk = getattribute(graph , :Zk)[end]
    n = getattribute(graph , :normalized)
    bound = n*lagrangeheuristic(graph)
    step = abs(Zk-bound)/(norm(res)^2)
    # First point
    α1 = 0
    z1 = Zk

    # Second point
    α2 = α
    z2 = αeval(α2,graph,bound)

    if (z2 - z1)/abs(z1) >= Δ
    # If second point gives an increase of more than Δ%, take it.
        return (λ += α2*step*res), bound
    elseif abs((z2 - z1)/z1) < Δ
    # If second point is similar to the first one, take the midpoint
        return (λ += α2/2*step*res), bound
    else
    # If second point is larger than first, take quarter-step
        if exhaustive
            return probingsubgradient(graph,λ,res,lagrangeheuristic,α/4,Δ;exhaustive=true)
        else
            return (λ += α2/4*step*res), bound
        end
    end
end

function cuttingplanes(graph,λ,res)
    ms = getattribute(graph , :lgmaster)
    Zk = getattribute(graph , :Zk)[end]
    nmult = getattribute(graph , :numlinks)

    λvar = getindex(ms, :λ)
    η = getindex(ms,:η)

    cut = @constraint(ms, η <= Zk + sum(λvar[j]*res[j] for j in 1:nmult))
    push!(getattribute(graph , :cuts), cut)

    solve(ms)
    return getvalue(λvar), getobjectivevalue(ms)
end


# Standard Lagrangean Heuristics
function fixbinaries(graph::ModelGraph,cat=[:Bin])
    if !haskey(graph.attributes,:mflat)
        setattribute(graph , :mflat, create_jump_graph_model(graph))
        getattribute(graph, :mflat).solver = getsolver(graph)
    end
    n = getattribute(graph , :normalized)
    mflat = getattribute(graph , :mflat)
    mflat.colVal = vcat([getmodel(n).colVal for n in getnodes(graph)]...)
    for j in 1:mflat.numCols
        if mflat.colCat[j] in cat
            mflat.colUpper[j] = mflat.colVal[j]
            mflat.colLower[j] = mflat.colVal[j]
        end
    end
    status = solve(mflat)
    if status == :Optimal
        return n*getobjectivevalue(mflat)
    else
        error("Heuristic model infeasible or unbounded")
    end
end

function fixintegers(graph::ModelGraph)
    fixbinaries(graph,[:Bin,:Int])
end
