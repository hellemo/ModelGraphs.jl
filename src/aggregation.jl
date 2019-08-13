#############################################################################################
# Aggregate Model
#############################################################################################
mutable struct AggregatedNode
    index::Int64
    obj_dict::Dict{Symbol,Any}
    variablemap::Dict{JuMP.VariableRef,JuMP.VariableRef}                    #map from aggregate model variable to original modelgraph variable
    constraintmap::Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}
    nl_constraintmap::Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}
    objective::Union{JuMP.AbstractJuMPScalar,Expr}                          #copy of original node objective
end
AggregatedNode(index::Int64) = AggregatedNode(index,Dict{Symbol,Any}(),Dict{JuMP.VariableRef,JuMP.VariableRef}(),
Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}(),Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}(),zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef}))

#Aggregation Info
mutable struct AggregationInfo
    nodes::Vector{AggregatedNode}
    linkvariables::Vector{VariableRef}
    linkconstraints::Vector{ConstraintRef}
    NLlinkconstraints::Vector{ConstraintRef}
end
AggregationInfo() = AggregationInfo(AggregatedNode[],VariableRef[],ConstraintRef[],ConstraintRef[])


#A JuMP model created from a aggregated ModelGraph
function AggregateModel()
    m = JuMP.Model()
    m.ext[:AggregationInfo] = AggregationInfo()
    return m
end

is_aggregate_model(m::JuMP.Model) = haskey(m.ext,:AggregationInfo) ? true : false  #check if the model is a graph model
assert_is_aggregate_model(m::JuMP.Model) = @assert is_aggregate_model(m)

getaggregationinfo(m::JuMP.Model) = haskey(m.ext, :AggregationInfo) ? m.ext[:AggregationInfo] : error("Model is not an aggregate model")
getlinkconstraints(m::JuMP.Model) = assert_aggregate_model(m) && getaggregationinfo(m).linkconstraints
getlinkvariables(m::JuMP.Model) = assert_aggregate_model(m) && getaggregationinfo(m).linkvariables
getNLlinkconstraints(m::JuMP.Model) = assert_aggregate_model(m) && getaggregationinfo(m).NLlinkconstraints
NHG.getnodes(m::JuMP.Model) = assert_aggregate_model(m) && getaggregationinfo(m).nodes

#Create a new new node on an AggregateModel
function add_aggregated_node!(m::JuMP.Model)
    assert_is_aggregate_model(m)
    i = getnumnodes(m)
    agg_node = AggregatedNode(i+1)
    push!(m.ext[:AggregationInfo].nodes,agg_node)
    return agg_node
end
getnumnodes(m::JuMP.Model) = length(getaggregationinfo(m).nodes)

JuMP.objective_function(node::AggregatedNode) = node.objective
JuMP.num_variables(node::AggregatedNode) = length(node.variablemap)
Base.getindex(node::AggregatedNode,s::Symbol) = node.obj_dict[s]
getaggnodevariables(node::AggregatedNode) = collect(keys(node.variablemap))
getaggnodeconstraints(node::AggregatedNode) = collect(keys(node.constraintmap))

#############################################################################################
# Aggregation Map
#############################################################################################
"""
    AggregationMap
    Mapping between variable and constraint reference of a node model to the Aggregated Model.
    The reference of the aggregated model can be obtained by indexing the map with the reference of the corresponding original modelnode.
"""
struct AggregationMap
    aggregate_model::JuMP.Model                             #An aggregate model
    varmap::Dict{JuMP.VariableRef,JuMP.VariableRef}         #map variables in original modelgraph to aggregatemodel
    conmap::Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}     #map constraints in original modelgraph to aggregatemodel
end

function Base.getindex(reference_map::AggregationMap, vref::JuMP.VariableRef)  #reference_map[node_var] --> aggregated_copy_var
    return reference_map.varmap[vref]
end

function Base.getindex(reference_map::AggregationMap, cref::JuMP.ConstraintRef)
    return reference_map.conmap[cref]
end
Base.broadcastable(reference_map::AggregationMap) = Ref(reference_map)

function Base.setindex!(reference_map::AggregationMap, graph_cref::JuMP.ConstraintRef,node_cref::JuMP.ConstraintRef)
    reference_map.conmap[node_cref] = graph_cref
end

function Base.setindex!(reference_map::AggregationMap, graph_vref::JuMP.VariableRef,node_vref::JuMP.VariableRef)
    reference_map.varmap[node_vref] = graph_vref
end

AggregationMap(m::JuMP.Model) = AggregationMap(m,Dict{JuMP.VariableRef,JuMP.VariableRef}(),Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}())

function Base.merge!(ref_map1::AggregationMap,ref_map2::AggregationMap)
    merge!(ref_map1.varmap,ref_map2.varmap)
    merge!(ref_map1.conmap,ref_map2.conmap)
end

#############################################################################################
# Aggregation Functions
#############################################################################################
#Aggregate modelgraph into AggregateModel
function aggregate(modelgraph::ModelGraph)
    aggregate_model = AggregateModel()
    reference_map = AggregationMap(aggregate_model)

    master_reference_map = _add_to_aggregate_model!(aggregate_model,getmastermodel(modelgraph))
    merge!(reference_map,master_reference_map)

    #COPY NODE MODELS INTO AGGREGATED MODEL
    has_nonlinear_objective = false                     #check if any nodes have nonlinear objectives
    for modelnode in getnodes(modelgraph)               #for each node in the model graph
        node_model = getmodel(modelnode)
        #Need to pass master reference so we use those variables instead of creating new ones
        node_reference_map = _add_to_aggregate_model!(aggregate_model,node_model)  #updates jump_graph_model,the jump_node, and the ref_map
        merge!(reference_map,node_reference_map)   #Update the reference_map

        #Check for nonlinear objective functions unless we know we already have one

        if has_nonlinear_objective != true
            has_nonlinear_objective = _has_nonlinear_obj(node_model)
        end
    end

    #OBJECTIVE FUNCTION
    if has_objective(modelgraph)
        agg_graph_obj = _copy_constraint_func(JuMP.objective_function(modelgraph),reference_map)
        JuMP.set_objective_function(aggregate_model,agg_graph_obj)
        JuMP.set_objective_sense(aggregate_model,JuMP.objective_sense(modelgraph))
    elseif has_NLobjective(modelgraph)
        #TODO
        # dgraph = JuMP.NLPEvaluator(modelgraph)
        # MOI.initialize(dgraph,[:ExprGraph])
        # graph_obj = MOI.objective_expr(dgraph)
        # _splice_nonlinear_variables!(graph_obj,reference_map)  #_splice_nonlinear_variables!(node_obj,var_maps[node])
        # JuMP.set_NL_objective(aggregate_model,JuMP.objective_sense(modelgraph,graph_obj))
    else
        set_node_objectives!(modelgraph,aggregate_model,reference_map,has_nonlinear_objective)          #set to the sum of the node objectives
    end

    #ADD LINK CONSTRAINTS
    for linkconstraint in all_linkconstraints(modelgraph)
        new_constraint = _copy_constraint(linkconstraint,reference_map)
        JuMP.add_constraint(aggregate_model,new_constraint)
    end

    #TODO ADD NLLINKCONSSTRAINTS
    # for nllinkconstraint in getallnllinkconstraints(modelgraph)
    # end

    return aggregate_model, reference_map
end

#Aggregate the subgraphs of a modelgrap where n_levels corresponds to how many levels remain, 0 means no subgraphs
function aggregate(graph::ModelGraph,n_levels::Int64)
    new_model_graph = ModelGraph()
    return new_model_graph
end


#previously _buildnodemodel!
function _add_to_aggregate_model!(aggregate_model::JuMP.Model,node_model::JuMP.Model)          #jump_node

    agg_node = add_aggregated_node!(aggregate_model)

    if JuMP.mode(node_model) == JuMP.DIRECT
        error("Cannot copy a node model in `DIRECT` mode. Use the `Model` ",
              "constructor instead of the `direct_model` constructor to be ",
              "able to aggregate into a new JuMP Model.")
    end

    #reference_map = GraphReferenceMap(m,MOIU.IndexMap())
    reference_map = AggregationMap(aggregate_model)

    #COPY VARIABLES
    for var in JuMP.all_variables(node_model)
        if is_linked_variable(var)                                       #if the variable is actually a link variable, we don't need to make a new one
            reference_map[var] = getlinkvariable(var)                    #get the master variable
        else
            new_x = JuMP.@variable(aggregate_model)                      #create an anonymous variable
            reference_map[var] = new_x                                   #map variable reference to new reference
            var_name = JuMP.name(var)
            new_name = var_name
            JuMP.set_name(new_x,new_name)
            if JuMP.start_value(var) != nothing
                JuMP.set_start_value(new_x,JuMP.start_value(var))
            end
            agg_node.variablemap[new_x] = var
        end
    end

    #COPY CONSTRAINTS
    #Use JuMP and check if I have a ScalarConstraint or VectorConstraint and use the reference map to create new constraints
    constraint_types = JuMP.list_of_constraint_types(node_model)
    for (func,set) in constraint_types
        constraint_refs = JuMP.all_constraints(node_model, func, set)
        for constraint_ref in constraint_refs
            constraint = JuMP.constraint_object(constraint_ref)
            new_constraint = _copy_constraint(constraint,reference_map)
            new_ref= JuMP.add_constraint(aggregate_model,new_constraint)
            agg_node.constraintmap[new_ref] = constraint_ref
            reference_map[constraint_ref] = new_ref
        end
    end

    #TODO Get nonlinear object data to work
    #COPY OBJECT DATA (JUMP CONTAINERS).  I don't really need this for this.  It would be nice for Aggregation though.
    for (name, value) in JuMP.object_dictionary(node_model)
        agg_node.obj_dict[name] = reference_map[value]
        #jump_node.obj_dict[name] = getindex.(reference_map, value)
    end

    #COPY NONLINEAR CONSTRAINTS
    nlp_initialized = false
    if node_model.nlp_data !== nothing
        d = JuMP.NLPEvaluator(node_model)           #Get the NLP evaluator object.  Initialize the expression graph
        MOI.initialize(d,[:ExprGraph])
        nlp_initialized = true
        for i = 1:length(node_model.nlp_data.nlconstr)
            expr = MOI.constraint_expr(d,i)                         #this returns a julia expression
            _splice_nonlinear_variables!(expr,node_model,reference_map)        #splice the variables from var_map into the expression
            new_nl_constraint = JuMP.add_NL_constraint(aggregate_model,expr)      #raw expression input for non-linear constraint
            constraint_ref = JuMP.ConstraintRef(node_model,JuMP.NonlinearConstraintIndex(i),new_nl_constraint.shape)
            #push!(jump_node.nl_constraints,constraint_ref)
            agg_node.nl_constraintmap[new_nl_constraint] = constraint_ref
            node_reference.nl_constraintmap[new_nl_constraint] = constraint_ref
            reference_map[constraint_ref] = new_nl_constraint
        end
    end

    #OBJECTIVE FUNCTION (store expression on aggregated_nodes)
    if !(_has_nonlinear_obj(node_model))
        #AFFINE OR QUADTRATIC OBJECTIVE
        new_objective = _copy_objective(node_model,reference_map)
        agg_node.objective = new_objective
    else
        #NONLINEAR OBJECTIVE
        if !nlp_initialized
            d = JuMP.NLPEvaluator(node_model)           #Get the NLP evaluator object.  Initialize the expression graph
            MOI.initialize(d,[:ExprGraph])
        end
        new_obj = _copy_nl_objective(d,variablemap)
        agg_node.objective = new_obj
    end
    return reference_map
end


#INTERNAL HELPER FUNCTIONS
function _set_node_objectives!(modelgraph::ModelGraph,aggregate_model::JuMP.Model,reference_map::AggregationMap,has_nonlinear_objective::Bool)
    if has_nonlinear_objective
        graph_obj = :(0) #NOTE Strategy: Build up a Julia expression (expr) and then call JuMP.set_NL_objective(expr)
        for node in getnodes(modelgraph)
            node_model = getmodel(node)
            JuMP.objective_sense(node_model) == MOI.MIN_SENSE ? sense = 1 : sense = -1
            d = JuMP.NLPEvaluator(node_model)
            MOI.initialize(d,[:ExprGraph])
            node_obj = MOI.objective_expr(d)
            _splice_nonlinear_variables!(node_obj,node_model,ref_map)  #_splice_nonlinear_variables!(node_obj,var_maps[node])
            node_obj = Expr(:call,:*,:($sense),node_obj)
            graph_obj = Expr(:call,:+,graph_obj,node_obj)  #update graph objective
        end
        JuMP.set_NL_objective(aggregate_model, MOI.MIN_SENSE, graph_obj)

    else
        graph_obj = sum(JuMP.objective_function(agg_node) for agg_node in getnodes(aggregate_model))    #NOTE: All of the node objectives are converted to Minimize (MOI.OptimizationSense(0))
        JuMP.set_objective(aggregate_model,MOI.MIN_SENSE,graph_obj)
    end
end

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
function _splice_nonlinear_variables!(expr::Expr,model::JuMP.Model,reference_map::AggregationMap)  #var_map needs to map the node_model index to the new model variable
    for i = 1:length(expr.args)
        if typeof(expr.args[i]) == Expr
            if expr.args[i].head != :ref             #keep calling _splice_nonlinear_variables! on the expression until it's a :ref. i.e. :(x[index])
                _splice_nonlinear_variables!(expr.args[i],model,reference_map)
            else  #it's a variable
                var_index = expr.args[i].args[2]     #this is the actual MOI index (e.g. x[1], x[2]) in the node model
                new_var = :($(reference_map.varmap[JuMP.VariableRef(model,var_index)]))
                expr.args[i] = new_var               #replace :(x[index]) with a :(JuMP.Variable)
            end
        end
    end
end

# COPY CONSTRAINT FUNCTIONS
function _copy_constraint_func(func::JuMP.GenericAffExpr,ref_map::AggregationMap)
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

function _copy_constraint_func(func::JuMP.GenericQuadExpr,ref_map::AggregationMap)
    new_aff = _copy_constraint_func(func.aff,ref_map)
    new_terms = OrderedDict([(JuMP.UnorderedPair(ref_map[pair.a],ref_map[pair.b]),coeff) for (pair,coeff) in func.terms])
    new_func = JuMP.GenericQuadExpr{Float64,JuMP.VariableRef}()
    new_func.terms = new_terms
    new_func.aff = new_aff
    #new_func.constant = func.constant
    return new_func
end

function _copy_constraint_func(func::JuMP.VariableRef,ref_map::AggregationMap)
    new_func = ref_map[func]
    return new_func
end

function _copy_constraint(constraint::JuMP.ScalarConstraint,ref_map::AggregationMap)
    new_func = _copy_constraint_func(constraint.func,ref_map)
    new_con = JuMP.ScalarConstraint(new_func,constraint.set)
    return new_con
end

function _copy_constraint(constraints::JuMP.VectorConstraint,ref_map::AggregationMap)
    new_funcs = [_copy_constraint_func(con.func) for con in constraints]
    new_con = JuMP.VectorConstraint(new_func,constraint.set,constraint.shape)
    return new_con
end

function _copy_constraint(constraint::LinkConstraint,ref_map::AggregationMap)
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
function _copy_objective(m::JuMP.Model,ref_map::AggregationMap)
    return _copy_objective(JuMP.objective_function(m),ref_map)
end

function _copy_objective(func::Union{JuMP.GenericAffExpr,JuMP.GenericQuadExpr},ref_map::AggregationMap)
    new_func = _copy_constraint_func(func,ref_map)
    return new_func
end

function _copy_objective(func::JuMP.VariableRef,ref_map::AggregationMap)
    new_func = ref_map[func]
    return new_func
end

function _copy_nl_objective(d::JuMP.NLPEvaluator,reference_map::AggregationMap)#variablemap::Dict{Int,VariableRef})
    new_obj = MOI.objective_expr(d)
    _splice_nonlinear_variables!(new_obj,d.m,reference_map)
    JuMP.objective_sense(d.m) == MOI.OptimizationSense(0) ? sense = 1 : sense = -1
    new_obj = Expr(:call,:*,:($sense),obj)
    return new_obj
end

function set_sum_of_objectives(graph::ModelGraph)
end


# #NOTE: This function can still be useful
# function _get_local_and_cross_links(model_graph::ModelGraph,nodes::Vector{ModelNode})
#     local_link_constraints = []  #all nodes in link constraint are in the set of nodes given
#     cross_link_constraints = []  #at least one node in a link constraint is not in the node subset (must be connected to the rest of the graph)
#     #Check each node's link constraints.  Add it to one of the above lists.
#     checked_links = []
#     for node in nodes
#         for (graph,links) in getlinkconstraints(node)  #NOTE:  Need to store linkconstraint references on nodes
#             for link in links
#                 if !(link in checked_links)  #If it's a new link constraint
#                     vars = link.terms.vars
#                     var_nodes = map(getnode,vars)
#                     if all(node -> node in nodes,var_nodes)
#                         push!(local_link_constraints,link)
#                     else
#                         push!(cross_link_constraints,link)
#                     end
#                     push!(checked_links,link)
#                 end
#             end
#         end
#     end
#     return local_link_constraints , cross_link_constraints
# end
#
# #TODO
# function get_local_and_shared_variables(model_graph::ModelGraph,edges::Vector{LinkingEdge})
# end


# function create_aggregate_graph(model_graph::ModelGraph,partition_data::PartitionData)
#     println("Building Aggregated Model Graph")
#
#     new_model_graph = ModelGraph()  #The new aggregated ModelGraph
#     partitions = partition_data.partitions
#     shared_entities = partition_data.shared_entities #Could be linkconstraints, shared variables, shared models, or pairs
#
#     #What kind of reference map do we need?  multiple?
#     variable_map = Dict{JuMP.VariableRef,JuMP.VariableRef}()
#     n_partitions = length(partitions)
#     #Aggregate Partitions
#     for i = 1:n_partitions  #create aggregate model for each partition
#         part = partitions[i]
#
#         #Get LinkingEdges from subgraphs if there are any.
#         local_shared_entities = partition_data.partition_entities[i]
#
#         aggregate_model,agg_ref_map = create_aggregate_model(model_graph,part,local_shared_entities)
#
#         #Update VariableMap
#         merge!(variable_map,agg_ref_map.varmap)
#
#         aggregate_node = add_node!(new_model_graph)
#         setmodel(aggregate_node,aggregate_model)
#
#     end
#
#     #NEW LINK CONSTRAINTS
#     for entity in shared_entities #Could be e.g. LinkConstraints
#         add_shared_entity!(new_model_graph,entity,variable_map)
#     end
#
#     return new_model_graph
#
# end
#
# function add_shared_entity!(graph::ModelGraph,edge::LinkingEdge,variable_map::Dict)
#     link_constraints = [LinkConstraint(link_ref) for link_ref in edge.linkconstraints]
#     for link_constraint in link_constraints
#         copy_constraint = AlgebraicGraphs._copy_constraint(link_constraint,variable_map)
#         JuMP.add_constraint(graph.linkmodel,copy_constraint)
#     end
# end
#
# function add_shared_entity!(graph::ModelGraph,graphvariable::JuMP.AbstractVariableRef,ref_map::GraphReferenceMap)
# end
#
# #Build up an aggregate model given a set of nodes.
# #function create_aggregate_model(model_graph::ModelGraph,nodes::Vector{ModelNode},link_edges::Vector{LinkingEdge})
# function create_aggregate_model(model_graph::ModelGraph,nodes::Vector,link_edges::Vector)
#     #local_links, cross_links = _get_local_and_cross_links(model_graph,nodes)
#
#     #Get corresponding ModelNodes and LinkConstraints for given indices
#     #Extract LinkConstraints from the Edges
#     link_constraints = []
#     for edge in link_edges
#         for link_ref in edge.linkconstraints
#             linkcon = LinkConstraint(link_ref)
#             push!(link_constraints,linkcon)
#         end
#     end
#
#     aggregate_model =  AlgebraicGraphs.JuMPGraphModel()         #Use a JuMPGraphModel so we can track the internal structure
#     jump_graph = getgraph(aggregate_model)
#
#     reference_map = GraphReferenceMap(aggregate_model,Dict(),Dict())
#     #variable_map = Dict()
#
#     has_nonlinear_objective = false
#
#     #COPY NODE MODELS INTO AGGREGATE MODEL
#     for model_node in nodes  #for each node in the model graph
#         nodeindex = getindex(model_graph,model_node)
#         jump_node = add_node!(aggregate_model,index = nodeindex)  #add at the same index
#
#         node_reference_map = _buildnodemodel!(aggregate_model,jump_node,model_node)
#
#         merge!(reference_map,node_reference_map)
#
#         node_model = getmodel(model_node)
#
#         if has_nonlinear_objective != true
#             has_nonlinear_objective = _has_nonlinear_obj(node_model)
#         end
#
#         #TODO Get nonlinear object data to work
#         #COPY OBJECT DATA (JUMP CONTAINERS).  I don't really need this for this.  It would be nice for Aggregation though.
#         for (name, value) in JuMP.object_dictionary(node_model)
#             #jump_node.obj_dict[name] = reference_map[value]
#             jump_node.obj_dict[name] = getindex.(reference_map, value)
#         end
#     end
#
#     #LOCAL LINK CONSTRAINTS
#     #TODO: Typing is causing issues here
#     for linkconstraint in link_constraints
#         new_constraint = _copy_constraint(linkconstraint,reference_map)
#         JuMP.add_constraint(aggregate_model,new_constraint)
#     end
#
#     #NODE OBJECTIVE
#     if !(has_nonlinear_objective)
#         graph_obj = sum(JuMP.objective_function(jump_node) for jump_node in getnodes(jump_graph)) #NOTE: All of the node objectives were converted to Minimize (MOI.OptimizationSense(0))
#         JuMP.set_objective(aggregate_model,MOI.OptimizationSense(0),graph_obj)
#     else
#         graph_obj = :(0) #NOTE Strategy: Build up a Julia expression (expr) and then call JuMP.set_NL_objective(expr)
#         for jump_node in getnodes(jump_graph)
#
#             #NOTE: Might be able to just convert AffExpr and QuadExpr into Julia Expressions to make this easier
#             id = getindex(jump_graph,jump_node)
#             node_model = getmodel(getnode(model_graph,id))
#             JuMP.objective_sense(node_model) == MOI.OptimizationSense(0) ? sense = 1 : sense = -1
#             d = JuMP.NLPEvaluator(node_model)
#             MOI.initialize(d,[:ExprGraph])
#             node_obj = MOI.objective_expr(d)
#             _splice_nonlinear_variables!(node_obj,ref_map)  #_splice_nonlinear_variables!(node_obj,var_maps[node])
#
#             node_obj = Expr(:call,:*,:($sense),node_obj)
#             graph_obj = Expr(:call,:+,graph_obj,node_obj)
#         end
#         JuMP.set_NL_objective(aggregate_model, MOI.OptimizationSense(0), graph_obj)
#     end
#
#     #Need to return multiple reference maps
#
#     return aggregate_model,reference_map
#
# end
