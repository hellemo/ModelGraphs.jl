JuMP.Model(modelgraph::ModelGraph;add_node_objectives = !(has_objective(model_graph))) = getmodel(aggregate(modelgraph,add_node_objectives = add_node_objectives))

function JuMP.optimize!(graph::ModelGraph,optimizer::JuMP.OptimizerFactory,kwargs...)
    println("Aggregating Model Graph...")
    aggregate_model,reference_map = aggregate(graph)
    println("Aggregation Complete")

    JuMP.optimize!(aggregate_model,optimizer;kwargs...)
    status = JuMP.termination_status(aggregate_model)

    # TODO Get all the correct status codes for copying a solution
    if JuMP.has_values(aggregate_model)
        _copysolution!(aggregate_model,graph)   #Now get our solution data back into the original ModelGraph
    end

    return status
end

function JuMP.optimize!(graph::ModelGraph,optimizer::AbstractModelGraphOptimizer,kwargs...)
    optimizer_model = optimizer.initialize_model(graph)
    status = optimize!(optimizer_model)
    _copysolution!(optimizer_model,graph)
    return status
end

# #TODO Make sure this still works. copy the solution from one graph to another where nodes and variables match
function _copysolution!(agg_model::JuMP.Model,model_graph::ModelGraph)
    for node in getnodes(agg_model)
        idx = getindex(jump_graph,node)

        model_node = getnode(model_graph,idx)
        #NOTE: This SHOULD work as long as Variable and Constraint Index always align
        for (jnodevar,modelvar) in node.variablemap
            model_node.variable_values[modelvar.index] = JuMP.value(jnodevar)
        end
        for (jnodeconstraint,modelconstraint) in node.constraintmap
            try
                model_node.constraint_dual_values[modelconstraint.index] = JuMP.dual(jnodeconstraint)
            catch ArgumentError #NOTE: Ipopt doesn't catch duals of quadtratic constraints
                continue
            end
        end
        for (jnodeconstraint,modelconstraint) in node.nl_constraintmap
            try
                model_node.nl_constraint_dual_values[modelconstraint.index] = JuMP.dual(jnodeconstraint)
            catch ArgumentError #NOTE: Ipopt doesn't catch duals of quadtratic constraints
                continue
            end
        end
    end

    #TODO Set dual solution values for LinkConstraints
end
