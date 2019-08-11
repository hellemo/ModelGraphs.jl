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
NestedHyperGraphs.getnodes(m::JuMP.Model) = assert_aggregate_model(m) && getaggregationinfo(m).nodes

function add_aggregated_node!(m::JuMP.Model)
    assert_is_aggregate_model(m)
    i = getnumnodes(m)
    agg_node = AggregatedNode(i+1)
    push!(m.ext[:AggregationInfo].nodes,agg_node)
    return agg_node
end


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
    aggregate_model::JuMP.Model   #An aggregate model
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
#function create_jump_graph_model(model_graph::AbstractModelGraph;add_node_objectives = !(hasobjective(model_graph)))  #Add objectives together if graph objective not provided
#Aggregate a graph into a node
function aggregate(graph::ModelGraph)#;levels = 0)  #levels = number of remaining subgraph levels
    aggmodel = AggregateModel()
    return aggmodel
end

#Aggregate the subgraphs of a modelgrap where n_levels corresponds to how many levels remain, 0 means no subgraphs
function aggregate(graph::ModelGraph,n_levels::Int64)
    new_model_graph = ModelGraph()
    return new_model_graph
end

#Aggregate a graph based on a model partition.  Return a new ModelGraph with possible subgraphs (If it was passed a recursive partition)
function aggregate(graph::ModelGraph,hyperpartition::HyperPartition)
    println("Building Aggregated Model Graph")

    new_model_graph = ModelGraph()

    #Get model subgraphs.  These will contain model nodes and LinkEdges.
    subgraphs_to_aggregate =  getsubgraphs(hyperpartitions.hypergraph_partitions)                     #map(x -> getnode(graph,x),hyperpartition.partition_nodes)  #partitions of modelnodes

    parent = getparent(hyperpartitions.hypergraph_partitions)

    shared_nodes = parent.sharednodes     #Could be linkconstraints, shared variables, shared models, or pairs
    shared_edges = parent.sharededges

    #What kind of reference map do we need?  multiple?
    variable_map = Dict{JuMP.VariableRef,JuMP.VariableRef}()

    #Aggregate "Leaf" Partitions
    for subgraph in subgraphs
        aggregate_model,agg_ref_map = aggregate(subgraph)
        merge!(variable_map,agg_ref_map.varmap)   #Update VariableMap
        aggregate_node = add_node!(new_model_graph)
        setmodel(aggregate_node,aggregate_model)
    end

    #NEW LINK CONSTRAINTS
    for shared_edge in shared_edges
        for linkconstraint in shared_edge.linkconstraints
            copy_constraint!(new_model_graph,linkconstraint,variable_map)
        end
    end

    #NEW LINK VARIABLES AND MASTER

    return new_model_graph

end



    return new_model_graph
end


#Aggregate modelgraph into AggregateModel
function aggregate(modelgraph::ModelGraph;add_node_objectives = true)
    aggregate_model = AggregateModel()
    reference_map = AggregationMap(aggregate_model)

    #TODO MASTER MODEL (LINK VARIABLES AND MASTER CONSTRAINTS)
    for linkvariable in getlinkvariables(modelgraph)
        #create new variables in agg_model
        master_reference_map = _add_to_aggregate_model!(aggregate_model,getmastermodel(modelgraph))
    end

    #COPY NODE MODELS INTO AGGREGATED MODEL
    has_nonlinear_objective = false                     #check if any nodes have nonlinear objectives
    for modelnode in getnodes(modelgraph)               #for each node in the model graph

        node_reference_map = _add_to_aggregate_model!(aggregate_model,model_node)  #updates jump_graph_model,the jump_node, and the ref_map
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
            graph_obj = sum(JuMP.objective_function(agg_node) for agg_node in getnodes(aggregate_model))    #NOTE: All of the node objectives were converted to Minimize (MOI.OptimizationSense(0))
            JuMP.set_objective(aggregate_model,MOI.MIN_SENSE,graph_obj)
        else
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
        end
    else #use the graph objecting
        graph_obj = objective_function(model_graph)
    end

    #ADD LINK CONSTRAINTS
    for linkconstraint in getalllinkconstraints(modelgraph)
        new_constraint = _copy_constraint(linkconstraint,reference_map)
        JuMP.add_constraint(aggregate_model,new_constraint)
    end

    return aggregate_model, reference_map
end

#previously _buildnodemodel!
function _add_to_aggregate_model!(aggregate_model::JuMP.Model,node_model::JuMP.Model)          #jump_node

    agg_node = add_agg_node!(m)

    if JuMP.mode(node_model) == JuMP.DIRECT
        error("Cannot copy a node model in `DIRECT` mode. Use the `Model` ",
              "constructor instead of the `direct_model` constructor to be ",
              "able to aggregate into a new JuMP Model.")
    end

    #reference_map = GraphReferenceMap(m,MOIU.IndexMap())
    reference_map = AggregationMap(aggregate_model)

    #COPY VARIABLES
    for var in JuMP.all_variables(node_model)

        # if isa(var,LinkVariableRef)
        #     #TODO: Check for LinkVariables
        # end

        new_x = JuMP.@variable(aggregate_model)    #create an anonymous variable
        reference_map[var] = new_x                                   #map variable reference to new reference
        var_name = JuMP.name(var)
        new_name = "$(getname(modelnode))"*var_name
        JuMP.set_name(new_x,new_name)
        if JuMP.start_value(var) != nothing
            JuMP.set_start_value(new_x,JuMP.start_value(var))
        end
        agg_node.variablemap[new_x] = var
    end

    #COPY CONSTRAINTS
    #Use JuMP and check if I have a ScalarConstraint or VectorConstraint and use the reference map to create new constraints
    #Using Option 1
    constraint_types = JuMP.list_of_constraint_types(node_model)
    for (func,set) in constraint_types
        constraint_refs = JuMP.all_constraints(node_model, func, set)
        for constraint_ref in constraint_refs
            constraint = JuMP.constraint_object(constraint_ref)
            new_constraint = _copy_constraint(constraint,reference_map)
            new_ref= JuMP.add_constraint(aggregate_model,new_constraint)
            agg_node.constraintmap[new_ref] = constraint_ref
            node_reference.constraintmap[new_ref] = constraint_ref
            reference_map[constraint_ref] = new_ref
        end
    end

    #TODO Get nonlinear object data to work
    #COPY OBJECT DATA (JUMP CONTAINERS).  I don't really need this for this.  It would be nice for Aggregation though.
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
            new_nl_constraint = JuMP.add_NL_constraint(aggregate_model,expr)      #raw expression input for non-linear constraint
            constraint_ref = JuMP.ConstraintRef(node_model,JuMP.NonlinearConstraintIndex(i),new_nl_constraint.shape)
            #push!(jump_node.nl_constraints,constraint_ref)
            agg_node.nl_constraintmap[new_nl_constraint] = constraint_ref
            node_reference.nl_constraintmap[new_nl_constraint] = constraint_ref
            reference_map[constraint_ref] = new_nl_constraint
        end
    end

    #OBJECTIVE FUNCTION
    if !(_has_nonlinear_obj(node_model))
        #AFFINE OR QUADTRATIC OBJECTIVE
        new_objective = _copy_objective(node_model,reference_map)
        agg_node.objective = new_objective
        node_reference.objective = new_objective
    else
        #NONLINEAR OBJECTIVE
        if !nlp_initialized
            d = JuMP.NLPEvaluator(node_model)           #Get the NLP evaluator object.  Initialize the expression graph
            MOI.initialize(d,[:ExprGraph])
        end
        new_obj = _copy_nl_objective(d,variablemap)
        agg_node.objective = new_obj
        node_reference.objective = new_obj
    end
    return reference_map
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
