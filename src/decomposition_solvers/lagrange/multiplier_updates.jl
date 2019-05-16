function subgradient!(lmodel::AbstractLagrangeModel,residual_equality::Vector,residual_inequality::Vector)
    #get an upper bound
    alpha = lmodel.solver_data.alpha
    lower_bound = lmodel.current_lower_bound

    if lmodel.solver_data.upper_bound_heurstic != nothing
        upper_bound = lmodel.solver_data.upper_bound_heurstic(lmodel)
        step = alpha*abs(lower_bound - upper_bound)/(norm([residual_equality;residual_inequality])^2)
    else
        upper_bound = nothing
        step = alpha
    end

    lambda_equality_delta = step*residual_equality  #update multipliers
    lambda_inequality_delta = step*residual_inequality

    return lambda_equality_delta,lambda_inequality_delta
end
