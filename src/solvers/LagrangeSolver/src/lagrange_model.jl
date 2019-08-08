abstract type AbstractLagrangeModel end

mutable struct LagrangeSolverData
    max_iterations::Int64
    epsilon::Float64
    time_limit::Float64
    alpha::Float64
    delta::Float64
    max_no_improvement::Int64

    function LagrangeSolverData()
        data = new()
        data.max_iterations = 100
        data.epsilon = 0.001
        data.time_limit = 3600
        data.alpha = 2.0
        data.delta = 0.5
        data.max_no_improvement = 10
        return data
    end
end

#IDEA: Create a LagrangeModel using ModelGraph information
mutable struct LagrangeModel <: AbstractLagrangeModel
    solver_data::LagrangeSolverData   #solver inputs
    solution::LagrangeSolution        #solution data

    #Callbacks we can expose
    multiplier_update_function::Function                        #function that takes a LagrangeModel argument and returns multiplier updates for equality and inequality multipliers
    upper_bound_function::Union{Nothing,Function}               #optional, can speed up convergence
    subproblem_solver::Union{Nothing,JuMP.OptimizerFactory}     #e.g. GLPK.Optimizer

    subproblems::Vector{JuMP.Model}                             #we assume these are all minimization problems

    equality_multipliers::Vector{Float64}                       #link equality constraint multipliers
    inequality_multipliers::Vector{Float64}                     #link inequality constraint multipliers

    link_eq_matrix::Union{Nothing,AbstractMatrix}               #link equality matrix
    link_ineq_matrix::Union{Nothing,AbstractMatrix}             #link inequality matrix

    b_eq::Vector{Float64}                                       #link equality right hand side
    b_ineq::Vector{Float64}                                     #link inequality right hand side

    link_eq_variables::Vector{JuMP.VariableRef}                #Vector of JuMP variables corresponding to link equality constraints
    link_ineq_variables::Vector{JuMP.VariableRef}               #Vector of JuMP variables corresponding to link inequality constraints

    current_lower_bound::Float64
    current_upper_bound::Float64

    residuals_equality::Vector           #A_eq*x_eq - b_eq
    residuals_inequality::Vector         #A_ineq*x_ineq - b_ineq
end

# NOTE: We can lift shared variables into equality constraints
# shared_variables::Dict{JuMP.AbstractJuMPScalar,JuMP.AbstractJuMPScalar} #these get dualized into simple linkconstraints
# shared_variable_multipliers::Vector{Float64}
#lagrange_model.shared_variables = getsharedvariables(graph)  #dictionary

LagrangeModel() = LagrangeModel(LagrangeSolverData(),LagrangeSolution(),subgradient!,nothing,nothing,
                                Vector{JuMP.Model}(),
                                Vector{Float64}(),Vector{Float64}(),nothing,nothing,Vector{Float64}(),Vector{Float64}(),
                                Vector{JuMP.VariableRef}(),Vector{JuMP.VariableRef}(),-Inf,Inf,
                                Vector{Float64}(),Vector{Float64}())

#Dual decomposition algorithm
#NOTE: This might get renamed to lagrange_solver depending how we deal with ADMM
function dual_decomposition_solve(graph::ModelGraph,args...;kwargs...)

    #NOTE: figure out which args to pass
    lagrange_model = LagrangeModel(graph)  #Lagrange model with data

    status = solve(lagrange_model)

    return status
end

#Populate a LagrangeModel using a ModelGraph
function LagrangeModel(graph::ModelGraph) #args,kwargs
    #Fill in all the data we need to do dual decomposition
    lagrange_model = LagrangeModel()

    #get data from graph
    lagrange_model.subproblems = [getmodel(node) for node in getnodes(graph)]
    link_eq_constraints = get_link_eq_constraints(graph)         #link equality constraints in graph
    link_ineq_constraints = get_link_ineq_constraints(graph)     #link inequality constraints in graph

    #Setup objective and multiplier list for each sub-problem
    for node in lagrange_model.subproblems
        prepare_subproblem!(node)
    end

    #Setup data structures
    link_eq_matrix,b_eq,link_eq_variables,link_eq_map = prepare_link_matrix(link_eq_constraints)
    link_ineq_matrix,b_ineq,link_ineq_variables,link_ineq_map = prepare_link_matrix(link_ineq_constraints)

    lagrange_model.link_eq_matrix = link_eq_matrix
    lagrange_model.b_eq = b_eq
    lagrange_model.link_eq_variables = link_eq_variables

    lagrange_model.link_ineq_matrix = link_ineq_matrix
    lagrange_model.b_ineq = b_ineq
    lagrange_model.link_ineq_variables = link_ineq_variables

    #Setup multipliers
    lagrange_model.equality_multipliers = zeros(size(link_eq_matrix,1))
    lagrange_model.inequality_multipliers = zeros(size(link_ineq_matrix,1))

    #TODO: Do the Multiplier warm start here?

    #Setup subproblem multiplier vectors
    _prepare_eq_multiplier_map!(link_eq_matrix,link_eq_variables)
    _prepare_ineq_multiplier_map!(link_ineq_matrix,link_ineq_variables)
    return lagrange_model
end

#Solve a Lagrange Model using dual decomposition or ADMM
function solve(lmodel::LagrangeModel)

    Zk_old = lmodel.current_lower_bound
    max_iterations = lmodel.solver_data.max_iterations
    epsilon = lmodel.solver_data.epsilon
    delta = lmodel.solver_data.delta

    no_improvement = 0
    for iter in 1:max_iterations
        println(iter)
        # Solve subproblems
        Zk = 0  #objective value
        for subproblem in lmodel.subproblems
           Zkn = solve_lagrange_subproblem!(lmodel,subproblem)
           Zk += Zkn                                      # add objective values
        end
        #steptaken = true

        lmodel.current_lower_bound = Zk  #Update lower bound
        # If no improvement in the lowerbound, increase the no improvement counter
        if Zk < Zk_old
            no_improvement += 1
        end
        Zk_old = Zk

        # If too many iterations have happened without improvement in the lower bound, decrease alpha (step-size)
        if no_improvement >= lmodel.solver_data.max_no_improvement
            no_improvement = 0
            lmodel.solver_data.alpha *= delta
        end

        A_eq = lmodel.link_eq_matrix
        A_ineq = lmodel.link_ineq_matrix

        x_eq = JuMP.value.(lmodel.link_eq_variables)
        x_ineq = JuMP.value.(lmodel.link_ineq_variables)

        # Update residuals for multplier calculation
        lmodel.residuals_equality = A_eq*x_eq - lmodel.b_eq

        #println(residuals_equality)
        lmodel.residuals_inequality = A_ineq*x_ineq - lmodel.b_ineq

        # Check convergence
        if norm(lmodel.residuals_equality) < epsilon && all(lmodel.residuals_inequality .>= 0)
            lmodel.solution.termination = :Optimal
            return :Optimal
        end

        #multiplier updates
        eq_multipliers_delta, ineq_multipliers_delta = lmodel.multiplier_update_function(lmodel)
        lmodel.equality_multipliers += eq_multipliers_delta
        lmodel.inequality_multipliers += ineq_multipliers_delta
        lmodel.inequality_multipliers[lmodel.inequality_multipliers .<= 0] .= 0  #NOTE: set to zero, or just don't update?
    end

    lmodel.solution.termination= :MaxIterationsReached
    return :MaxIterationsReached

end

#Use current multiplier values to solve a lagrange sub-problem
function solve_lagrange_subproblem!(lagrange_model::LagrangeModel,subproblem::JuMP.Model)
    # Restore original objective function.  We need to remake it to update multipliers
    #IDEA: I could use parameter JuMP and make the multipliers into parameters

    m = subproblem

    eq_multipliers = lagrange_model.equality_multipliers
    ineq_multipliers = lagrange_model.inequality_multipliers

    scale = m.ext[:objective_scale]
    JuMP.set_objective_function(m,scale*m.ext[:original_objective])

    #Add dualized part to objective function for each linkconstraint that hits this node
    obj = JuMP.objective_function(m)

    #Add equality multiplier terms
    for term in m.ext[:eq_multipliers]  #could have duplicate terms
        coeff = term[1]
        var = term[2]
        multiplier_index = term[3]
        obj += coeff*var*eq_multipliers[multiplier_index]
    end

    #Add inequality multiplier terms
    for term in m.ext[:ineq_multipliers]
        coeff = term[1]
        var = term[2]
        multiplier_index = term[3]
        obj += coeff*var*ineq_multipliers[multiplier_index]
    end

    #node[:lagrange_objective] = obj
    JuMP.set_objective_function(m,obj)

    JuMP.optimize!(m)

    m.ext[:current_objective_value] = JuMP.objective_value(m)
    return JuMP.objective_value(m)
end

function _prepare_eq_multiplier_map!(link_matrix::SparseMatrixCSC{Float64,Int64},link_variables::Vector)
    for i in 1:size(link_matrix)[1] #loop through rows
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
    for i in 1:size(link_matrix)[1] #loop through rows
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


function admm_solve(graph::ModelGraph)
end
# data.cp_upper_bound = 1e6
# data.cutting_plane_solver = nothing
# data.subproblem_solver = nothing

# data.upper_bound_function = nothing
# data.multiplier_update_method = subgradient!
