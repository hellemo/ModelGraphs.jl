mutable struct LagrangeSolverData
    max_iterations::Int
    multiplier_update_method::Function
    epsilon::Float64
    time_limit::Float64
    alpha::Float64
    heuristic::Function  #e.g. fix binaries, fix master problem solution
    delta::Float64
    max_no_improvement::Int64
    cp_upper_bound::Float64  #cutting plane upperbound if using cutting plane update method
    subproblem_solver::Union{Nothing,JuMP.OptimizerFactory}  #If nothing, use the optimizer set on the subproblem
    cutting_plane_solver::Union{Nothing,JuMP.OptimizerFactory}

    function LagrangeSolverData()
        data = new()
        data.max_iterations = 10
        data.update_method = SubGradient()
        data.epsilon = 0.001
        data.timelimit = 3600
        data.alpha = 2
        data.heurstic = nothing
        data.delta = 0.5
        data.cutting_plane_solver = nothing
        data.node_solver = nothing  #ClpSolver()
        return data
    end
end

#IDEA: Convert ModelGraph into LagrangeModel
#NOTE: Need to handle possibility of interval linkconstraints
mutable struct LagrangeModel
    solver_data::LagrangeSolverData
    solution::LagrangeSolution

    subproblems::Vector{JuMP.Model}   #we assume they are all minimization problems

    shared_variables::Dict{JuMP.AbstractJuMPScalar,JuMP.AbstractJuMPScalar} #these get dualized into simple linkconstraints

    link_eq_constraints::Vector{JuMP.AbstractConstraint}   #Ax = b
    link_ineq_constraints::Vector{JuMP.AbstractConstraint} #Ax <= b

    shared_variable_multipliers::Vector{Float64}

    equality_multipliers::Vector{Float64}

    inequality_multipliers::Vector{Float64}

    # link_eq_matrix::AbstractMatrix
    # link_ineq_matrix::AbstractMatrix
    # eq_link_variables::Vector{JuMP.AbstractJuMPScalar}
    # ineq_link_variables::Vector{JuMP.AbstractJuMPScalar}

end


LagrangeModel() = LagrangeModel(LagrangeSolverData(),LagrangeSolution(),Vector{JuMP.Model},
                                Dict{JuMP.AbstractJuMPScalar,JuMP.AbstractJuMPScalar}(),Vector{JuMP.AbstractConstraint}(),Vector{JuMP.AbstractConstraint}(),
                                Vector{Float64}(),Vector{Float64}(),Vector{Float64}())

function LagrangeModel(graph::ModelGraph)
    #Fill in all the data we need to do dual decomposition
    subproblems = [getmodel(node) for node in getnodes(graph)]
    shared_variables = getsharedvariables(graph)  #dictionary


end

#Flip objective, add dictionary entry for multipliers for this model
function prepare_subproblem!(m::JuMP.Model)
    m.ext[:multiplier_map] = Dict()
    m.ext[:objective_scale] = 1
    m = getmodel(node)
    obj = JuMP.objective_function(m)
    if JuMP.objective_sense(m) == MOI.MAX_SENSE
        m.ext[:objective_scale] = -1
        JuMP.set_objective_sense(m,MOI.MIN_SENSE)
    end
    m.ext[:original_objective] = obj
end


#Dual decomposition algorithm
function dual_decomposition_solve(graph::ModelGraph,args...;kwargs...)

    # #Initialize multipliers

    #TODO return mapping that maps variable terms to link constraints (multipliers)
    link_eq_matrix,link_eq_variables,b_eq = prepare_link_eq_matrix(graph)
    link_ineq_matrix,link_ineq_variables,b_ineq = prepare_link_ineq_matrix(graph)

    #NOTE: Initialize multiplier vectors here

    lagrange_model = LagrangeModel(graph)  #Lagrange model with data

    #Setup LagrangeModel
    equality_multipliers = zeros(size(link_eq_matrix,1))
    inequality_multipliers = zeros(size(link_ineq_matrix,2))

    # if initialize_multipliers == :lp_relaxation
    #     solve_lp_relaxation!(graph)
    #     #TODO: Update lagrange model multipliers
    # end

    #Setup multiplier map for each sub-problem
    for node in lagrange_model.subproblems
        prepare_subproblem!(node)
    end

    #Setup multiplier maps
    #TODO create a proper link constraint mapping
    for link_con in getlinkconstraints(graph)
        for (var,coeff) in link.func.terms
            node = getnode(var)
            term = (var,coeff)
            node[:multiplier_map][term] = index  #Need to get reference to correct linkconstraint here
        end
    end

    starttime = time()

    #Solve Lagrange Model
    ### ITERATIONS ###
    for iter in 1:max_iterations

        # Solve subproblems
        Zk = 0  #objective value
        for subprob in lagrange_model.subproblems
           Zkn = solve_lagrange_subproblem!(node)
           Zk += Zkn                                      # add objective values
        end
        steptaken = true

        # If no improvement in the lowerbound, increase the no improvement counter
        if Zk < Zk_old
            no_improvement += 1
        end
        Zk_old = Zk

        # If too many iterations have happened without improvement in the lower bound, decrease alpha
        if no_improvement >= maximum_no_improvement
            no_improvement = 0
            alpha *= delta
        end

        # Update residuals for multplier calculation
        residuals_equality = Pi_eq*value.(link_eq_variables) - b_eq
        residuals_inequality = Pi_ineq*value.(link_ineq_variables) - b_ineq

        # Check convergence
        if norm(residuals_equality) < epsilon && norm(residuals_inequality) < epsilon
            return :Optimal
        end

        # # Check time limit
        # if tstamp > timelimit
        #     return :TimeLimit
        # end

        # Update multipliers alpha,Zk,lambda,residual,lagrangeheuristic
        multipliers, upper_bound = subgradient(alpha,Zk,multipliers,residuals_equality,heur)
    end
end

function admm_solve(graph::ModelGraph)
end
