
JuMP.Model(modelgraph::ModelGraph;add_node_objectives = !(has_objective(model_graph))) = getmodel(aggregate(modelgraph,add_node_objectives = add_node_objectives))



function jump_solve(graph::ModelGraph,optimizer::JuMP.OptimizerFactory;scale = 1.0,kwargs...)
    println("Building JuMP Model...")
    m_jump,reference_map = create_jump_graph_model(graph)
    println("Finished Creating JuMP Model")

    #Reset the scaled objective
    #TODO Get rid of this with an actual graph objective
    JuMP.set_objective_function(m_jump,scale*JuMP.objective_function(m_jump))

    JuMP.optimize!(m_jump,optimizer;kwargs...)
    status = JuMP.termination_status(m_jump)

    graph.jump_model = m_jump

    # TODO Get all the correct status codes for copying a solution
    if JuMP.has_values(m_jump)
        _copysolution!(getgraph(m_jump),graph)          #Now get our solution data back into the original ModelGraph
    end

    return status
end

#TODO Remove scale argument.  Allow direct interface with the graph objective function
function JuMP.optimize!(graph::AbstractModelGraph,optimizer::JuMP.OptimizerFactory;scale = 1.0,kwargs...)
    status = jump_solve(graph,optimizer,scale = scale,kwargs...)
    return status
end

# #TODO Make sure this still works. copy the solution from one graph to another where nodes and variables match
function _copysolution!(agg_model::JuMP.Model,model_graph::ModelGraph)
    for node in getnodes(jump_graph)
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
