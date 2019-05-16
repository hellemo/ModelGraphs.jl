function subgradient!(lmodel::AbstractLagrangeModel)
    #get an upper bound
    alpha = lmodel.solver_data.alpha
    lower_bound = lmodel.current_lower_bound

    if lmodel.upper_bound_function != nothing
        lmodel.upper_bound = lmodel.upper_bound_function(lmodel)
        step = alpha*abs(lower_bound - lmodel.upper_bound)/(norm([lmodel.residual_equality;lmodel.residual_inequality])^2)
    else
        upper_bound = nothing
        step = alpha
    end

    lambda_equality_delta = step*lmodel.residuals_equality  #update multipliers
    lambda_inequality_delta = step*lmodel.residuals_inequality

    return lambda_equality_delta,lambda_inequality_delta
end

function cutting_plane!(lmodel::AbstractLagrangeModel)
end
