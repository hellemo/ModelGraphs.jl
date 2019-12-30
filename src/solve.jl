JuMP.Model(modelgraph::ModelGraph;add_node_objectives = !(has_objective(model_graph))) = getmodel(aggregate(modelgraph,add_node_objectives = add_node_objectives))

function JuMP.optimize!(graph::ModelGraph,optimizer::JuMP.OptimizerFactory;kwargs...)
    println("Aggregating Model Graph...")
    aggregate_model,reference_map = aggregate(graph)

    println("Optimizing Aggregated Model...")
    JuMP.optimize!(aggregate_model,optimizer;kwargs...)
    status = JuMP.termination_status(aggregate_model)

    if JuMP.has_values(aggregate_model)     # TODO Get all the correct status codes for copying a solution
        _copysolution!(graph,reference_map)   #Now get our solution data back into the original ModelGraph
        println("Found Solution")
    end

    return status
end

# function JuMP.optimize!(graph::ModelGraph,optimizer::AbstractModelGraphOptimizer,kwargs...)
#     optimizer_model = initialize_model(optimizer,graph)
#     status = optimize!(optimizer_model)
#     _copysolution!(optimizer_model,graph)
#     return status
# end

# #TODO Make sure this still works. copy the solution from one graph to another where nodes and variables match
function _copysolution!(modelgraph::ModelGraph,ref_map::AggregationMap)
    for node in all_nodes(modelgraph)
        for var in JuMP.all_variables(node)
            node.variable_values[var] = JuMP.value(ref_map[var])
        end
    end

    #TODO Copy constraint duals
    # for (jnodeconstraint,modelconstraint) in node.constraintmap
    #     try
    #         model_node.constraint_dual_values[modelconstraint.index] = JuMP.dual(jnodeconstraint)
    #     catch ArgumentError #NOTE: Ipopt doesn't catch duals of quadtratic constraints
    #         continue
    #     end
    # end
    #     for (jnodeconstraint,modelconstraint) in node.nl_constraintmap
    #         try
    #             model_node.nl_constraint_dual_values[modelconstraint.index] = JuMP.dual(jnodeconstraint)
    #         catch ArgumentError #NOTE: Ipopt doesn't catch duals of quadtratic constraints
    #             continue
    #         end
    #     end
    # end

    # #TODO Copy linkconstraint duals

end
