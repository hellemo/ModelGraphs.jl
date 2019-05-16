mutable struct LagrangeSolverData
    max_iterations::Int
    multiplier_update_method::Function
    epsilon::Float64
    time_limit::Float64
    alpha::Float64
    delta::Float64
    max_no_improvement::Int64

    cp_upper_bound::Float64  #cutting plane upperbound if using cutting plane update method
    subproblem_solver::Union{Nothing,JuMP.OptimizerFactory}  #If nothing, use the optimizer set on the subproblem
    cutting_plane_solver::Union{Nothing,JuMP.OptimizerFactory}

    upperbound_heuristic::Function  #e.g. fix binaries, fix master problem solution
    multiplier_update::Function

    function LagrangeSolverData()
        data = new()
        data.max_iterations = 10
        data.update_method = subgradient
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

    # link_eq_constraints::Vector{JuMP.AbstractConstraint}   #Ax = b
    # link_ineq_constraints::Vector{JuMP.AbstractConstraint} #Ax <= b

    #shared_variable_multipliers::Vector{Float64}


    equality_multipliers::Vector{Float64}
    inequality_multipliers::Vector{Float64}

    link_eq_matrix::AbstractMatrix
    b_eq::Vector{Float64}

    link_ineq_matrix::AbstractMatrix
    b_ineq::Vector{Float64}

    link_eq_variables::Vector{JuMP.VariableRef}
    link_ineq_variables::Vector{JuMP.VariableRef}

end


LagrangeModel() = LagrangeModel(LagrangeSolverData(),LagrangeSolution(),Vector{JuMP.Model},
                                Dict{JuMP.AbstractJuMPScalar,JuMP.AbstractJuMPScalar}(),Vector{JuMP.AbstractConstraint}(),Vector{JuMP.AbstractConstraint}(),
                                Vector{Float64}(),Vector{Float64}(),Vector{Float64}())

#Convert ModelGraph to a LagrangeModel
function LagrangeModel(graph::ModelGraph)
    #Fill in all the data we need to do dual decomposition
    lagrange_model = LagrangeModel()

    lagrange_model.subproblems = [getmodel(node) for node in getnodes(graph)]

    link_eq_constraints = get_link_eq_constraints(graph)         #equality
    link_ineq_constraints = get_link_ineq_constraints(graph)

    #Setup objective and multiplier map for each sub-problem
    for node in lagrange_model.subproblems
        prepare_subproblem!(node)
    end

    #lagrange_model.shared_variables = getsharedvariables(graph)  #dictionary

    #multiplier_map = Dict() #linkconstraint --> multiplier index

    #Setup link data structures
    link_eq_matrix,b_eq,link_eq_variables,link_eq_map = prepare_link_matrix(link_eq_constraints)
    link_ineq_matrix,b_ineq,link_ineq_variables,link_ineq_map = prepare_link_matrix(link_ineq_constraints)

    lagrange_model.link_eq_matrix = link_eq_matrix
    lagrange_model.b_eq = b_eq
    lagrange_model.link_eq_variables = link_eq_variables

    lagrange_model.link_ineq_matrix = link_eq_matrix
    lagrange_model.b_ineq = b_ineq
    lagrange_model.link_ineq_variables = link_ineq_variables

    #Setup LagrangeModel
    lagrange_model.equality_multipliers = zeros(size(link_eq_matrix,1))
    lagrange_model.inequality_multipliers = zeros(size(link_ineq_matrix,2))

    #TODO: Do the Multiplier warm start here

    _prepare_eq_multiplier_map!(link_eq_matrix,link_eq_variables)
    _prepare_ineq_multiplier_map!(link_ineq_matrix,link_ineq_variables)
    #setup subproblem multiplier maps

end

function _prepare_eq_multiplier_map!(link_matrix::SparseMatrixCSC{Float64,Int64},link_variables::Vector)
    for i in 1:size(link_eq_matrix)[1] #loop through rows
        row = link_matrix[i,:]
        for j in row.nzind
            coeff = link_matrix[i,j]
            var = link_variables[j]
            model = var.model
            push!(model.ext[:eq_multipliers],(coeff,var,i))
        end
    end
end

function _prepare_ineq_multiplier_map!(link_matrix::SparseMatrixCSC{Float64,Int64},link_variables::Vector)
    for i in 1:size(link_eq_matrix)[1] #loop through rows
        row = link_matrix[i,:]
        for j in row.nzind
            coeff = link_matrix[i,j]
            var = link_variables[j]
            model = var.model
            push!(model.ext[:ineq_multipliers],(coeff,var,i))
        end
    end
end


#Flip objective, add dictionary entry for multipliers for this model
function prepare_subproblem!(m::JuMP.Model)
    m.ext[:eq_multipliers] = []
    m.ext[:ineq_multipliers] = []
    m.ext[:objective_scale] = 1
    obj = JuMP.objective_function(m)
    if JuMP.objective_sense(m) == MOI.MAX_SENSE
        m.ext[:objective_scale] = -1
        JuMP.set_objective_sense(m,MOI.MIN_SENSE)
    end
    m.ext[:original_objective] = obj
end

#Solve a Lagrange Model using dual decomposition or ADMM
function solve(model::LagrangeModel)
    start_time = time()

    for iter in 1:max_iterations
        # Solve subproblems
        Zk = 0  #objective value
        for subproblem in lagrange_model.subproblems
           Zkn = solve_lagrange_subproblem!(subproblem)
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


        Pi_eq = model.link_eq_matrix
        Pi_ineq = model.link_ineq_matrix
        # Update residuals for multplier calculation
        residuals_equality = Pi_eq*value.(model.link_eq_variables) - model.b_eq
        residuals_inequality = Pi_ineq*value.(model.link_ineq_variables) - model.b_ineq

        # Check convergence
        if norm(residuals_equality) + norm(residuals_inequality) < epsilon
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

#Use current multiplier values to solve a lagrange sub-problem
function solve_lagrange_subproblem!(lagrange_model::LagrangeModel,subproblem::JuMP.Model)
    # Restore original objective function.  We need to remake it to update multipliers
    #IDEA: I could use parameter JuMP and make the multipliers into parameters
    scale = m.ext[:objective_scale]
    JuMP.set_objective_function(m,scale*m.ext[:original_objective])
    # Add dualized part to objective function for each linkconstraint that hits this node
    multiplier_map = m.ext[:multiplier_map]  #mapping: lambda --> variable and its coefficient

    #Setup Dual Decomposition objective
    obj = JuMP.objective_function(m)

    #Add equality multiplier terms
    for (term,multiplier_index) in m.ext[:eq_multiplier_map]  #could have duplicate terms
        coeff = term[2]
        var = term[1]
        obj += coeff*var*eq_multipliers[multiplier_index]
    end

    #Add inequality multiplier terms
    for (term,multiplier_index) in m.ext[:ineq_multiplier_map]
        coeff = term[2]
        var = term[1]
        obj += coeff*var*ineq_multipliers[multiplier_index]
    end

    node[:lagrange_objective] = obj
    JuMP.set_objective_function(m,obj)

    JuMP.optimize!(m)

    m.ext[:current_objective_value] = JuMP.objective_value(m)
    return JuMP.objective_value(m)
end


#Dual decomposition algorithm
function dual_decomposition_solve(graph::ModelGraph,args...;kwargs...)

    # #Initialize multipliers
    #NOTE: Initialize multiplier vectors here
    lagrange_model = LagrangeModel(graph)  #Lagrange model with data
    #TODO return mapping that maps variable terms to link constraints (multipliers)
    if initialize_multipliers == :lp_relaxation
        solve_lp_relaxation!(graph)
        #TODO: Update lagrange model multipliers
    end

    status = solve(lagrange_model)

end

function admm_solve(graph::ModelGraph)
end
