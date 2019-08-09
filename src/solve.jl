"""
create_jump_graph_model
    Create a JuMP model through an aggregation procedure on the ModelNodes, LinkConstraints, GraphVariables, and GraphConstraints.
    Return a JuMPGraphModel which is a Model with graph extension data containing a JuMPGraph.
"""
#THE IDEA : Create a JuMP model by aggregating all of the nodes in a ModelGraph together

JuMP.Model(modelgraph::ModelGraph;add_node_objectives = !(hasobjective(model_graph))) = getmodel(aggregate(modelgraph,add_node_objectives = add_node_objectives))

create_jump_graph_model(modelgraph::ModelGraph) = aggregate(modelgraph)


#function create_jump_graph_model(model_graph::AbstractModelGraph;add_node_objectives = !(hasobjective(model_graph)))  #Add objectives together if graph objective not provided
#IDEA: Aggregate graph into a single node
function aggregate(model_graph::ModelGraph)
    # jump_graph_model = JuMPGraphModel()
    # jump_graph = StructureGraphs.copy_graph_to(model_graph,to_graph_type = JuMPGraph)  #Create empty JuMP graph with nodes and edges
    # jump_graph_model.ext[:Graph] = jump_graph
    aggregate_model = AggregationModel()
    reference_map = AggregationMap(aggregate_model)

    # COPY NODE MODELS INTO AGGREGATED MODEL
    has_nonlinear_objective = false                     #check if any nodes have nonlinear objectives
    for model_node in getnodes(model_graph)             #for each node in the model graph
        nodeindex = getindex(model_graph,model_node)

        jump_node = getnode(jump_graph,nodeindex)

        node_reference_map = _buildnodemodel!(jump_graph_model,jump_node,model_node)  #updates jump_graph_model,the jump_node, and the ref_map

        #Update the reference_map
        merge!(reference_map,node_reference_map)

        #Check for nonlinear objective functions unless we know we already have one
        node_model = getmodel(model_node)
        if has_nonlinear_objective != true
            has_nonlinear_objective = _has_nonlinear_obj(node_model)
        end
    end

    #OBJECTIVE FUNCTION
    if add_node_objectives
        if !(has_nonlinear_objective)
            graph_obj = sum(JuMP.objective_function(jump_node) for jump_node in getnodes(jump_graph)) #NOTE: All of the node objectives were converted to Minimize (MOI.OptimizationSense(0))
            JuMP.set_objective(jump_graph_model,MOI.OptimizationSense(0),graph_obj)
        else
            graph_obj = :(0) #NOTE Strategy: Build up a Julia expression (expr) and then call JuMP.set_NL_objective(expr)
            for jump_node in getnodes(jump_graph)

                #NOTE: Might be able to just convert AffExpr and QuadExpr into Julia Expressions to make this easier
                id = getindex(jump_graph,jump_node)
                node_model = getmodel(getnode(model_graph,id))
                JuMP.objective_sense(node_model) == MOI.OptimizationSense(0) ? sense = 1 : sense = -1
                d = JuMP.NLPEvaluator(node_model)
                MOI.initialize(d,[:ExprGraph])
                node_obj = MOI.objective_expr(d)
                _splice_nonlinear_variables!(node_obj,node_model,ref_map)  #_splice_nonlinear_variables!(node_obj,var_maps[node])

                node_obj = Expr(:call,:*,:($sense),node_obj)
                graph_obj = Expr(:call,:+,graph_obj,node_obj)
            end
            JuMP.set_NL_objective(jump_graph_model, MOI.OptimizationSense(0), graph_obj)
        end
    else
        #TODO. Use GraphObjective
    end

    #LINK CONSTRAINTS
    for linkconstraint in get_all_linkconstraints(model_graph)
        new_constraint = _copy_constraint(linkconstraint,reference_map)
        JuMP.add_constraint(jump_graph_model,new_constraint)
    end

    #TODO Link VARIABLES  #Should do this before copying model nodes
    for graph_variable in getgraphvariables(model_graph)
    end

    #TODO Master Model
    for graph_constraint in getgraphconstraints(model_graph)
    end

    return jump_graph_model, reference_map
end

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
function _copysolution!(jump_graph::JuMPGraph,model_graph::ModelGraph)
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

    #TODO Set dual solution values for Graph Constraints (LinkConstraints)




end

#INTERNAL HELPER FUNCTIONS
function _has_nonlinear_obj(m::JuMP.Model)
    if m.nlp_data != nothing
        if m.nlp_data.nlobj != nothing
            return true
        end
    end
    return false
end

#COPY OBJECTIVE FUNCTIONS
#splice variables into a constraint expression
function _splice_nonlinear_variables!(expr::Expr,model::JuMP.Model,reference_map::GraphReferenceMap)  #var_map needs to map the node_model index to the new model variable
    for i = 1:length(expr.args)
        if typeof(expr.args[i]) == Expr
            if expr.args[i].head != :ref   #keep calling _splice_nonlinear_variables! on the expression until it's a :ref. i.e. :(x[index])
                _splice_nonlinear_variables!(expr.args[i],model,reference_map)
            else  #it's a variable
                var_index = expr.args[i].args[2]     #this is the actual MOI index (e.g. x[1], x[2]) in the node model
                #new_var = :($(reference_map.index_map.varmap[var_index].value))  #Get MOIIndex.value
                #Need index --> new variable
                new_var = :($(reference_map.varmap[JuMP.VariableRef(model,var_index)]))
                #new_var = :($(var_map[var_index]))   #get the JuMP variable from var_map using the index
                expr.args[i] = new_var               #replace :(x[index]) with a :(JuMP.Variable)
            end
        end
    end
end

# COPY CONSTRAINT FUNCTIONS
function _copy_constraint_func(func::JuMP.GenericAffExpr,ref_map::GraphReferenceMap)
    terms = func.terms
    new_terms = OrderedDict([(ref_map[var_ref],coeff) for (var_ref,coeff) in terms])
    new_func = JuMP.GenericAffExpr{Float64,JuMP.VariableRef}()
    new_func.terms = new_terms
    new_func.constant = func.constant
    return new_func
end

function _copy_constraint_func(func::JuMP.GenericAffExpr,var_map::Dict{JuMP.VariableRef,JuMP.VariableRef})
    terms = func.terms
    new_terms = OrderedDict([(var_map[var_ref],coeff) for (var_ref,coeff) in terms])
    new_func = JuMP.GenericAffExpr{Float64,JuMP.VariableRef}()
    new_func.terms = new_terms
    new_func.constant = func.constant
    return new_func
end

function _copy_constraint_func(func::JuMP.GenericQuadExpr,ref_map::GraphReferenceMap)
    new_aff = _copy_constraint_func(func.aff,ref_map)
    new_terms = OrderedDict([(JuMP.UnorderedPair(ref_map[pair.a],ref_map[pair.b]),coeff) for (pair,coeff) in func.terms])
    new_func = JuMP.GenericQuadExpr{Float64,JuMP.VariableRef}()
    new_func.terms = new_terms
    new_func.aff = new_aff
    #new_func.constant = func.constant
    return new_func
end

function _copy_constraint_func(func::JuMP.VariableRef,ref_map::GraphReferenceMap)
    new_func = ref_map[func]
    return new_func
end

function _copy_constraint(constraint::JuMP.ScalarConstraint,ref_map::GraphReferenceMap)
    new_func = _copy_constraint_func(constraint.func,ref_map)
    new_con = JuMP.ScalarConstraint(new_func,constraint.set)
    return new_con
end

function _copy_constraint(constraints::JuMP.VectorConstraint,ref_map::GraphReferenceMap)
    new_funcs = [_copy_constraint_func(con.func) for con in constraints]
    new_con = JuMP.VectorConstraint(new_func,constraint.set,constraint.shape)
    return new_con
end

function _copy_constraint(constraint::LinkConstraint,ref_map::GraphReferenceMap)
    new_func = _copy_constraint_func(constraint.func,ref_map)
    new_con = JuMP.ScalarConstraint(new_func,constraint.set)
    return new_con
end

function _copy_constraint(constraint::LinkConstraint,var_map::Dict{JuMP.VariableRef,JuMP.VariableRef})
    new_func = _copy_constraint_func(constraint.func,var_map)
    new_con = JuMP.ScalarConstraint(new_func,constraint.set)
    return new_con
end

#COPY OBJECTIVE FUNCTIONS
function _copy_objective(m::JuMP.Model,ref_map::GraphReferenceMap)
    return _copy_objective(JuMP.objective_function(m),ref_map)
end

function _copy_objective(func::Union{JuMP.GenericAffExpr,JuMP.GenericQuadExpr},ref_map::GraphReferenceMap)
    new_func = _copy_constraint_func(func,ref_map)
    return new_func
end

function _copy_objective(func::JuMP.VariableRef,ref_map::GraphReferenceMap)
    new_func = ref_map[func]
    return new_func
end

function _copy_nl_objective(d::JuMP.NLPEvaluator,reference_map::GraphReferenceMap)#variablemap::Dict{Int,VariableRef})
    new_obj = MOI.objective_expr(d)
    _splice_nonlinear_variables!(new_obj,d.m,reference_map)
    JuMP.objective_sense(d.m) == MOI.OptimizationSense(0) ? sense = 1 : sense = -1
    new_obj = Expr(:call,:*,:($sense),obj)
    return new_obj
end
