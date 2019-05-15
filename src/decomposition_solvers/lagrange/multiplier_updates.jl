function subgradient(alpha,Zk,lambda,residual,lagrange_heuristic)
    #get an upper bound
    if lagrange_heuristic() != nothing
        upper_bound = lagrange_heuristic()
        step = alpha*abs(Zk-upper_bound)/(norm(residual)^2)
    else
        upper_bound = nothing
        step = alpha
    end
    lambda += step*residual  #update multipliers
    return lambda,upper_bound
end
