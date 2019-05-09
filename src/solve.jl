#TODO IndexMap with NonlinearConstraintIndex

"""
GraphReferenceMap
    Mapping between variable and constraint reference of a node model to the Aggregated Model.
    The reference of the aggregated model can be obtained by indexing the map with the reference of the corresponding reference of the original node model.
"""
#TODO Make sure this is working correctly
struct GraphReferenceMap
    model::JuMP.Model   #An aggregate model
    varmap::Dict{JuMP.VariableRef,JuMP.VariableRef}
    conmap::Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}
end
function Base.getindex(reference_map::GraphReferenceMap, vref::JuMP.VariableRef)  #reference_map[node_var] --> aggregated_copy_var
    #return JuMP.VariableRef(reference_map.model,reference_map.index_map[JuMP.index(vref)])
    return reference_map.varmap[vref]
end
function Base.getindex(reference_map::GraphReferenceMap, cref::JuMP.ConstraintRef)
    #return JuMP.ConstraintRef(reference_map.model,reference_map.index_map[JuMP.index(cref)],cref.shape)
    return reference_map.conmap[cref]
end
Base.broadcastable(reference_map::GraphReferenceMap) = Ref(reference_map)

function Base.setindex!(reference_map::GraphReferenceMap, graph_cref::JuMP.ConstraintRef,node_cref::JuMP.ConstraintRef)
    reference_map.conmap[node_cref] = graph_cref
    #reference_map.index_map.conmap[node_cref.index] = graph_cref.index
end
function Base.setindex!(reference_map::GraphReferenceMap, graph_vref::JuMP.VariableRef,node_vref::JuMP.VariableRef)
    reference_map.varmap[node_vref] = graph_vref
    #reference_map.index_map.varmap[node_vref.index] = graph_vref.index
end

GraphReferenceMap(m::JuMP.Model) = GraphReferenceMap(m,Dict{JuMP.VariableRef,JuMP.VariableRef}(),Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}())

function Base.merge!(ref_map1::GraphReferenceMap,ref_map2::GraphReferenceMap)
    merge!(ref_map1.varmap,ref_map2.varmap)
    merge!(ref_map1.conmap,ref_map2.conmap)
    # for (k,v) in ref_map2.index_map.varmap
    #     ref_map1.index_map.varmap[k] = v
    # end
    # for (k,v) in ref_map2.index_map.conmap
    #     ref_map1.index_map.conmap[k] = v
    # end
end

"""
create_jump_graph_model
    Create a JuMP model through an aggregation procedure on the ModelNodes, LinkConstraints, GraphVariables, and GraphConstraints.
    Return a JuMPGraphModel which is a Model with graph extension data containing a JuMPGraph.
"""
#THE IDEA : Create a JuMP model by aggregating all of the nodes in a ModelGraph together
function create_jump_graph_model(model_graph::AbstractModelGraph;add_node_objectives = !(hasobjective(model_graph)))  #Add objectives together if graph objective not provided
    println("Creating JuMP Model")
    jump_graph_model = JuMPGraphModel()

    jump_graph = StructureGraphs.copy_graph_to(model_graph,to_graph_type = JuMPGraph)  #Create empty JuMP graph with nodes and edges
    jump_graph_model.ext[:Graph] = jump_graph

    reference_map = GraphReferenceMap(jump_graph_model)

    # COPY NODE MODELS INTO AGGREGATED MODEL
    println("Aggregating Models")
    has_nonlinear_objective = false                     #check if any nodes have nonlinear objectives
    for model_node in getnodes(model_graph)             #for each node in the model graph
        nodeindex = getindex(model_graph,model_node)
        #println(nodeindex)
        jump_node = getnode(jump_graph,nodeindex)

        node_reference_map = _buildnodemodel!(jump_graph_model,jump_node,model_node)  #updates jump_graph_model,the jump_node, and the ref_map

        #Update the reference_map
        #NOTE: this might be too much data to store
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

    println("Adding Link Constraints")
    #LINK CONSTRAINTS
    for linkconstraint in get_all_linkconstraints(model_graph)
        new_constraint = _copy_constraint(linkconstraint,reference_map)
        JuMP.add_constraint(jump_graph_model,new_constraint)
    end

    #TODO GRAPH VARIABLES
    for graph_variable in getgraphvariables(model_graph)
    end

    #TODO GRAPH CONSTRAINTS
    for graph_constraint in getgraphconstraints(model_graph)
    end

    return jump_graph_model, reference_map
end

#Function to build a node model for a flat graph model
function _buildnodemodel!(m::JuMP.Model,jump_node::JuMPNode,model_node::ModelNode)
    node_model = getmodel(model_node)

    if JuMP.mode(node_model) == JuMP.DIRECT
        error("Cannot copy a node model in `DIRECT` mode. Use the `Model` ",
              "constructor instead of the `direct_model` constructor to be ",
              "able to aggregate into a new JuMP Model.")
    end

    #reference_map = GraphReferenceMap(m,MOIU.IndexMap())
    reference_map = GraphReferenceMap(m)
    #COPY VARIABLES
    for var in JuMP.all_variables(node_model)
        new_x = JuMP.@variable(m)    #create an anonymous variable

        reference_map[var] = new_x                                   #map variable reference to new reference
        var_name = JuMP.name(var)
        new_name = "$(getlabel(jump_node))$(getindex(getgraph(m),jump_node))."*var_name
        JuMP.set_name(new_x,new_name)
        if JuMP.start_value(var) != nothing
            JuMP.set_start_value(new_x,JuMP.start_value(var))
        end
        push!(jump_node.variablelist,new_x)
    end

    #COPY CONSTRAINTS
    #NOTE Option 1: Currently implemented
    #Use JuMP and check if I have ScalarConstraint or VectorConstraint and use my index map to create new constraints

    #Option 2: Not implemented
    #Use MOI to copy which would hit all of these constraints.  See MOI.copy_to for ideas about how to do this.

    #Using Option 1
    constraint_types = JuMP.list_of_constraint_types(node_model)
    for (func,set) in constraint_types
        constraint_refs = JuMP.all_constraints(node_model, func, set)
        for constraint_ref in constraint_refs
            constraint = JuMP.constraint_object(constraint_ref)
            new_constraint = _copy_constraint(constraint,reference_map)

            #reference_map[constraint] = new_constraint

            ref= JuMP.add_constraint(m,new_constraint)
            push!(jump_node.constraintlist,ref)
        end
    end

    #TODO Get nonlinear to work
    #COPY OBJECT DATA (JUMP CONTAINERS)
    # for (name, value) in JuMP.object_dictionary(node_model)
    #     jump_node.obj_dict[name] = reference_map[value]
    #     #jump_node.obj_dict[name] = getindex.(reference_map, value)
    # end

    #COPY NONLINEAR CONSTRAINTS
    nlp_initialized = false
    if node_model.nlp_data !== nothing
        d = JuMP.NLPEvaluator(node_model)           #Get the NLP evaluator object.  Initialize the expression graph
        MOI.initialize(d,[:ExprGraph])
        nlp_initialized = true
        for i = 1:length(node_model.nlp_data.nlconstr)
            expr = MOI.constraint_expr(d,i)                         #this returns a julia expression
            _splice_nonlinear_variables!(expr,node_model,reference_map)        #splice the variables from var_map into the expression
            #println(expr)
            new_nl_constraint = JuMP.add_NL_constraint(m,expr)      #raw expression input for non-linear constraint
            constraint_ref = JuMP.ConstraintRef(node_model,JuMP.NonlinearConstraintIndex(i),new_nl_constraint.shape)
            #reference_map[constraint_ref] = new_nl_constraint
            push!(jump_node.constraintlist,new_nl_constraint)
        end
    end

    #OBJECTIVE FUNCTION
    if !(_has_nonlinear_obj(node_model))
        #AFFINE OR QUADTRATIC OBJECTIVE
        new_objective = _copy_objective(node_model,reference_map)
        jump_node.objective = new_objective
    else
        #NONLINEAR OBJECTIVE
        if !nlp_initialized
            d = JuMP.NLPEvaluator(node_model)           #Get the NLP evaluator object.  Initialize the expression graph
            MOI.initialize(d,[:ExprGraph])
        end
        new_obj = _copy_nl_objective(d,variablemap)
        jump_node.objective = new_obj
    end
    return reference_map
end

#Create a JuMP model and solve with a MPB compliant solver
#buildjumpmodel!(graph::AbstractModelGraph) = graph.serial_model = create_jump_graph_model(graph)

function jump_solve(graph::ModelGraph,optimizer::JuMP.OptimizerFactory;scale = 1.0,kwargs...)
    println("Aggregating Models...")
    m_jump,reference_map = create_jump_graph_model(graph)
    println("Finished Creating JuMP Model")

    #Reset the scaled objective
    #TODO Get rid of this with an actual graph objective
    JuMP.set_objective_function(m_jump,scale*JuMP.objective_function(m_jump))

    JuMP.optimize!(m_jump,optimizer;kwargs...)
    status = JuMP.termination_status(m_jump)

    #TODO Get all the correct status codes for copying a solution
    # if JuMP.has_values(m_jump)
    #     _copysolution!(getgraph(m_jump),graph,reference_map)          #Now get our solution data back into the original ModelGraph
    # end
    graph.jump_model = m_jump
    return status
end

#check if graph has a MOI connected solver
#TODO Remove scale argument.  Allow direct interface with the graph objective function
function JuMP.optimize!(graph::AbstractModelGraph,optimizer::JuMP.OptimizerFactory;scale = 1.0,kwargs...)
    status = jump_solve(graph,optimizer,scale = scale,kwargs...)
    return status
end

# #TODO Make sure this still works. copy the solution from one graph to another where nodes and variables match
# function _copysolution!(jump_graph::JuMPGraph,model_graph::ModelGraph,ref_map::GraphReferenceMap)
#     for node in getnodes(jump_graph) #graph 1 is the JuMPGraph
#         idx = getindex(jump_graph,node)
#         model_node = getnode(model_graph,idx)       #get the corresponding node or edge in graph2
#
#         for variable in node.variablelist
#             model_node_var = ref_map[variable]
#
#             #JuMP.set_start_value(model_node_var,JuMP.value(variable))
#         end
#
#     end
#     #TODO Set dual solution values for constraints
# end

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
#TODO Fix this
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
#function _copy_constraint_func(func::JuMP.GenericAffExpr{Float64,V <: JuMP.AbstractVariableRef},ref_map::GraphReferenceMap)
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
    new_aff = copy_constraint_func(func.aff)
    new_terms = OrderedDict([(ref_map[var_ref],coeff) for (var_ref,coeff) in terms])
    new_func = JuMP.GenericQuadExpr{Float64,JuMP.VariableRef}()
    new_func.terms = new_terms
    new_func.aff = new_aff
    new_func.constant = func.constant
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

# function _add_to_objective!(m::JuMP.Model,add_obj::JuMP.AbstractJuMPScalar)
#     obj = JuMP.objective_function(m)
#     new_obj = _add_to_objective!(m,obj,add_obj)
#     JuMP.set_objecive_function(m,new_obj)
#     return new_obj
# end
#
# function _add_to_objective!(m::JuMP.Model,obj::JuMP.AbstractJuMPScalar,add_obj::JuMP.AbstractJuMPScalar)
#     new_obj = obj + add_obj
#     new_obj::JuMP.AbstractJuMPScalar
#     return new_obj
# end
#
# function _add_to_objective!(m::JuMP.Model,obj::JuMP.AbstractJuMPScalar,expr::Expr)
#     JuMP.set_NL_objective(jump_model, :Min, obj)
# end

# #IDEA: Use MOI to create a backend copy and then fill in the missing JuMP information
# function my_copy_to(aggregate_model::JuMP.Model, src_model::JuMP.Model, copy_names::Bool)
#
#     dest = backend(aggregate_model)
#     src = backend(src_model)
#
#     #MOI.empty!(dest)
#
#     idxmap = MOI.Utilities.IndexMap()
#
#     #Copy variables
#     vis_src = MOI.get(src, MOI.ListOfVariableIndices())
#     vars = MOI.add_variables(dest, length(vis_src))
#     for (vi, var) in zip(vis_src, vars)
#         idxmap.varmap[vi] = var
#     end
#
#     # Copy variable attributes
#     MOI.Utilities.pass_attributes(dest, src, copy_names, idxmap, vis_src)
#
#     # Copy model attributes
#     MOI.Utilities.pass_attributes(dest, src, copy_names, idxmap)
#
#     # Copy constraints
#     for (F, S) in MOI.get(src, MOI.ListOfConstraints())
#         # do the rest in copyconstraints! which is type stable
#         MOI.Utilities.copyconstraints!(dest, src, copy_names, idxmap, F, S)
#     end
#
#     return idxmap
# end
