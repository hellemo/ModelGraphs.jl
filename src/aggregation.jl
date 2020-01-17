#############################################################################################
# Aggregate: IDEA: Group nodes together into a larger node
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
getlinkconstraints(m::JuMP.Model) = is_aggregate_model(m) && getaggregationinfo(m).linkconstraints
getlinkvariables(m::JuMP.Model) = is_aggregate_model(m) && getaggregationinfo(m).linkvariables
getNLlinkconstraints(m::JuMP.Model) = is_aggregate_model(m) && getaggregationinfo(m).NLlinkconstraints
getnodes(m::JuMP.Model) = is_aggregate_model(m) && getaggregationinfo(m).nodes

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
    Mapping between variable and constraint reference of a ModelGraph to an Aggregated Model.
    The reference of the aggregated model can be obtained by indexing the map with the reference of the corresponding original modelnode.
"""
struct AggregationMap
    aggregate_model::JuMP.AbstractModel                             #An aggregate model (Could be another ModelGraph)
    varmap::Dict{JuMP.VariableRef,JuMP.VariableRef}                 #map variables in original modelgraph to aggregatemodel
    conmap::Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}             #map constraints in original modelgraph to aggregatemodel
end

function Base.getindex(reference_map::AggregationMap, vref::JuMP.VariableRef)  #reference_map[node_var] --> aggregated_copy_var
    return reference_map.varmap[vref]
end

function Base.getindex(reference_map::AggregationMap, lvref::LinkVariableRef)  #reference_map[node_var] --> aggregated_copy_var
    return reference_map.varmap[lvref.vref]
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

# function Base.setindex!(reference_map::AggregationMap, graph_vref::JuMP.VariableRef,node_vref::JuMP.VariableRef)
#     reference_map.varmap[node_vref] = graph_vref
# end

AggregationMap(m::JuMP.AbstractModel) = AggregationMap(m,Dict{JuMP.VariableRef,JuMP.VariableRef}(),Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}())

function Base.merge!(ref_map1::AggregationMap,ref_map2::AggregationMap)
    merge!(ref_map1.varmap,ref_map2.varmap)
    merge!(ref_map1.conmap,ref_map2.conmap)
end

#############################################################################################
# Aggregation Functions
#############################################################################################
#Aggregate modelgraph into AggregateModel
#Group ModelGraph to a ModelNode
convert_to_node(modelgraph::ModelGraph) = aggregate(modelgraph)


function aggregate(modelgraph::ModelGraph)  #group, collapse,
    aggregate_model = AggregateModel()
    reference_map = AggregationMap(aggregate_model)

    #TODO: Get rid of master node idea
    master_node = getmasternode(modelgraph)
    _add_to_aggregate_model!(aggregate_model,getmodel(master_node),reference_map)

    #COPY NODE MODELS INTO AGGREGATED MODEL
    has_nonlinear_objective = false                      #check if any nodes have nonlinear objectives
    for modelnode in all_nodes(modelgraph)               #for each node in the model graph
        node_model = getmodel(modelnode)

        #Need to pass master reference so we use those variables instead of creating new ones
        _add_to_aggregate_model!(aggregate_model,node_model,reference_map)  #updates aggregate_model and reference_map

        #Check for nonlinear objective functions unless we know we already have one
        if has_nonlinear_objective != true
            has_nonlinear_objective = _has_nonlinear_obj(node_model)
        end
    end

    #OBJECTIVE FUNCTION
    if !(has_objective(modelgraph)) && !has_nonlinear_objective
        _set_node_objectives!(modelgraph)  #set modelgraph objective function
        _set_node_objectives!(modelgraph,aggregate_model,reference_map,has_nonlinear_objective) #set aggregate_model objective function
    end

    if has_objective(modelgraph)
        agg_graph_obj = _copy_constraint_func(JuMP.objective_function(modelgraph),reference_map)
        JuMP.set_objective_function(aggregate_model,agg_graph_obj)
        JuMP.set_objective_sense(aggregate_model,JuMP.objective_sense(modelgraph))
    # elseif has_NLobjective(modelgraph)
    #     #TODO
    #     error("NL graph objective not yet supported on a ModelGraph")
    #     # dgraph = JuMP.NLPEvaluator(modelgraph)
    #     # MOI.initialize(dgraph,[:ExprGraph])
    #     # graph_obj = MOI.objective_expr(dgraph)
    #     # _splice_nonlinear_variables!(graph_obj,reference_map)  #_splice_nonlinear_variables!(node_obj,var_maps[node])
    #     # JuMP.set_NL_objective(aggregate_model,JuMP.objective_sense(modelgraph,graph_obj))
    # else
    #     _set_node_objectives!(modelgraph,aggregate_model,reference_map,has_nonlinear_objective)  #Set objective on the aggregate model
    end

    #ADD LINK CONSTRAINTS
    for linkconstraint in all_linkconstraints(modelgraph)
        new_constraint = _copy_constraint(linkconstraint,reference_map)
        JuMP.add_constraint(aggregate_model,new_constraint)
    end

    #TODO ADD NLLINKCONSSTRAINTS
    # for nllinkconstraint in getallnllinkconstraints(modelgraph)
    # end

    modelnode = ModelNode()
    set_model(modelnode,aggregate_model)

    # return aggregate_model, reference_map
    return modelnode,reference_map
end

#Aggregate a graph using a model partition.  Return a new ModelGraph with possible subgraphs (If it was passed a recursive partition), link constraints and link variables
#IDEA: Create new ModelGraph with subgraphs based on partition object.
#Group subgraphs together for solver interface
function aggregate(graph::ModelGraph,hyperpartition::Partition,hypermap::Dict)
    println("Creating Partitioned ModelGraph...")

    #Create New ModelGraphs
    parent_dict = Dict()
    for parent in hyperpartition.partitionroots
        new_model_graph = ModelGraph()
        parent_dict[parent] = new_model_graph
    end

    top_model_graph = parent_dict[hyperpartition.partitionroots[1]]
    reference_map = AggregationMap(top_model_graph)  #old model graph => new modelgraph


    #ADD BOTTOM LEVEL NODES: Aggregate subgraphs and then create bottom level nodes
    #submodelgraphs = []
    for partition in hyperpartition.leafpartitions

        #hypergraph = partition.hypergraph
        hypernodes = partition.hypernodes
        hyperedges = partition.hyperedges
        modelnodes = ModelNode[hypermap[node] for node in hypernodes]
        linkedges = LinkEdge[hypermap[edge] for edge in hyperedges]

        submodelgraph = ModelGraphs.induced_modelgraph(modelnodes,linkedges)
        #push!(submodelgraphs,submodelgraph)

        aggregate_node,agg_ref_map = aggregate(submodelgraph) #creates new model

        merge!(reference_map,agg_ref_map)

        parent_graph = parent_dict[partition.parent]

        #aggregate_node = add_node!(parent_graph)
        #set_model(aggregate_node,aggregate_model)
        add_node!(parent_graph,aggregate_node)
    end

    #Now add shared nodes and shared edges
    for parent in hyperpartition.partitionroots
        shared_nodes = parent.sharednodes     #shared nodes (modelnodes)
        shared_edges = parent.sharededges     #shared edges (linkedges)

        parent_mg = parent_dict[parent]

        # LINK VARIABLES (MASTER NODE)
        # master = aggregate(shared_nodes) #get linkvariables from shared nodes
        # set_master(parent_mg,master)
        #master = Model()
        for shared_node in shared_nodes
            error("Shared nodes not supported yet")
            #identify edges here and figure out which link variables to make
        end
        #parent_mg.masternode = master

        #LINK CONSTRAINTS
        for shared_edge in shared_edges
            #linkedge = findlinkedge(graph,shared_edge)
            linkedge = hypermap[shared_edge]
            for linkconstraintref in linkedge.linkconstraints
                linkconstraint = LinkConstraint(linkconstraintref)
                new_con = _copy_constraint(linkconstraint,reference_map)
                JuMP.add_constraint(parent_mg,new_con)  #add linkconstraint to the modelgraph
            end
        end

        if parent.parent != nothing
            parent_subgraph = parent_dict[parent.parent]
            add_subgraph!(subgraph,new_model_graph)
        end
    end
    return top_model_graph,reference_map
end

#TODO
#Aggregate the subgraphs of a modelgraph where n_levels corresponds to how many levels remain, 0 means no subgraphs
# function aggregate(graph::ModelGraph,n_levels::Int64)
#     new_model_graph = ModelGraph()
#     return new_model_graph
# end

#Create a modelgraph where the subgraphs are grouped into ModelNodes
function group_subgraphs(graph::ModelGraph)
end

function _add_to_aggregate_model!(aggregate_model::JuMP.Model,node_model::JuMP.Model,aggregation_map::AggregationMap)

    agg_node = add_aggregated_node!(aggregate_model)

    if JuMP.mode(node_model) == JuMP.DIRECT
        error("Cannot copy a node model in `DIRECT` mode. Use the `Model` ",
              "constructor instead of the `direct_model` constructor to be ",
              "able to aggregate into a new JuMP Model.")
    end

    reference_map = AggregationMap(aggregate_model)
    constraint_types = JuMP.list_of_constraint_types(node_model)
    #COPY VARIABLES
    for var in JuMP.all_variables(node_model)
        if is_linked_variable(var)                                       #if the variable is actually a link variable, we don't need to make a new one
            reference_map[var] = aggregation_map[getlinkvariable(var)]   #get the master variable
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

    #COPY ALL OTHER CONSTRAINTS
    #Use JuMP and check if I have a ScalarConstraint or VectorConstraint and use the reference map to create new constraints
    for (func,set) in constraint_types
        constraint_refs = JuMP.all_constraints(node_model, func, set)
        for constraint_ref in constraint_refs
            constraint = JuMP.constraint_object(constraint_ref)
            if func == JuMP.VariableRef
                var = constraint.func
                if is_linked_variable(var)  #Don't add multiple copies of a linked variable constraint
                    continue
                end
            end
            new_constraint = _copy_constraint(constraint,reference_map)
            new_ref= JuMP.add_constraint(aggregate_model,new_constraint)
            agg_node.constraintmap[new_ref] = constraint_ref
            reference_map[constraint_ref] = new_ref
        end

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
            agg_node.nl_constraintmap[new_nl_constraint] = constraint_ref
            reference_map[constraint_ref] = new_nl_constraint
        end
    end

    #TODO Get nonlinear object data to work
    #COPY OBJECT DATA (JUMP CONTAINERS).  I don't really need this for this.  It would be nice for Aggregation though.
    for (name, value) in JuMP.object_dictionary(node_model)
        #agg_node.obj_dict[name] = reference_map[value]
        if typeof(value) in [JuMP.VariableRef,JuMP.ConstraintRef,LinkVariableRef]
            agg_node.obj_dict[name] = getindex.(reference_map, value)
        end
    end

    #OBJECTIVE FUNCTION (store expression on aggregated_nodes)
    if !(_has_nonlinear_obj(node_model))
        #AFFINE OR QUADTRATIC OBJECTIVE
        new_objective = _copy_objective(node_model,reference_map)
        #TODO: Check objective sense
        agg_node.objective = new_objective
    else
        #NONLINEAR OBJECTIVE
        if !nlp_initialized
            d = JuMP.NLPEvaluator(node_model)           #Get the NLP evaluator object.  Initialize the expression graph
            MOI.initialize(d,[:ExprGraph])
        end
        new_obj = _copy_nl_objective(d,reference_map)
        agg_node.objective = new_obj
    end

    merge!(aggregation_map,reference_map)

    return reference_map
end

#Set aggregate model objective to sum of ModelGraph node objectives
function _set_node_objectives!(modelgraph::ModelGraph,aggregate_model::JuMP.Model,reference_map::AggregationMap,has_nonlinear_objective::Bool)
    if has_nonlinear_objective
        graph_obj = :(0) #NOTE Strategy: Build up a Julia expression (expr) and then call JuMP.set_NL_objective(expr)
        for node in all_nodes(modelgraph)
            node_model = getmodel(node)
            JuMP.objective_sense(node_model) == MOI.MIN_SENSE ? sense = 1 : sense = -1
            d = JuMP.NLPEvaluator(node_model)
            MOI.initialize(d,[:ExprGraph])
            node_obj = MOI.objective_expr(d)
            _splice_nonlinear_variables!(node_obj,node_model,reference_map)  #_splice_nonlinear_variables!(node_obj,var_maps[node])
            node_obj = Expr(:call,:*,:($sense),node_obj)
            graph_obj = Expr(:call,:+,graph_obj,node_obj)  #update graph objective
        end
        JuMP.set_NL_objective(aggregate_model, MOI.MIN_SENSE, graph_obj)
    else
        #TODO: Fix issue with setting maximize
        graph_obj = sum(JuMP.objective_function(agg_node) for agg_node in getnodes(aggregate_model))    #NOTE: All of the node objectives are converted to Minimize (MOI.OptimizationSense(0))
        JuMP.set_objective(aggregate_model,MOI.MIN_SENSE,graph_obj)
    end
end

function _set_node_objectives!(modelgraph::ModelGraph)
    # graph_obj = zero(JuMP.GenericAffExpr{Float64, JuMP.VariableRef})
    graph_obj = zero(JuMP.GenericQuadExpr{Float64, JuMP.VariableRef})  #testing changing this to quadratic expression
    for node in all_nodes(modelgraph)
        sense = JuMP.objective_sense(node)
        s = sense == MOI.MAX_SENSE ? -1.0 : 1.0
        JuMP.add_to_expression!(graph_obj,s,JuMP.objective_function(node))
    end

    JuMP.set_objective(modelgraph,MOI.MIN_SENSE,graph_obj)
end

#Create a ModelGraph from a set of ModelNodes and LinkEdges
#IDEA: Make copies instead?
function induced_modelgraph(modelnodes::Vector{ModelNode},linkedges::Vector{LinkEdge})
    submg = ModelGraph()
    for node in modelnodes
        push!(submg.modelnodes,node)
        submg.node_idx_map[node] = length(submg.modelnodes)
    end
    link_idx = 0
    for linkedge in linkedges
        #TODO: Make sure linkedge nodes actually connect the modelnodes
        push!(submg.linkedges,linkedge)
        submg.edge_idx_map[linkedge] = length(submg.linkedges)
        submg.linkedge_map[linkedge.nodes] = linkedge
        for linkconstraintref in linkedge.linkconstraints
            link_idx += 1
            # idx = linkconstraintref.idx #these can be duplicates with subgraphs
            linkconstraint = LinkConstraint(linkconstraintref)
            submg.linkconstraints[link_idx] = linkconstraint
        end
    end
    return submg
end

# Create an induced subgraph from a set of nodes
# function induced_modelgraph(modelgraph::ModelGraph,modelnodes::Vector{ModelNode})
#     #figure out supporting edges using hypergraph
# end

# #Create a ModelGraph from a given Hypergraph
# function create_sub_modelgraph(modelgraph::ModelGraph,hypergraph::HyperGraph,hyper_map::Dict)
#     submg = ModelGraph()
#
#     for hypernode in getnodes(hypergraph)
#         modelnode = hyper_map[hypernode]
#         push!(submg.modelnodes,modelnode)
#         submg.node_idx_map[modelnode] = length(graph.modelnodes)
#         #modelnode = getnode(modelgraph,hypernode)
#         #submg.modelnodes[hypernode] = modelnode
#     end
#
#     i = 1
#     for hyperedge in gethyperedges(hypergraph)
#         #linkedge = findlinkedge(modelgraph,hyperedge)  #this could be in a subgraph
#         linkedge = hyper_map[hyperedge]
#         add_edge!(submg,linkedge)
#         #submg.linkedges[hyperedge] = linkedge
#         for linkconstraintref in linkedge.linkconstraints
#             linkconstraint = LinkConstraint(linkconstraintref)
#             submg.linkconstraints[i] = linkconstraint
#             i += 1
#         end
#     end
#     return submg
# end

# #ADD VARIABLE CONSTRAINTS SEPARATELY.  Need to do this because linked variables get counted multiple times.
# for (func,set) in constraint_types
#     if func == JuMP.VariableRef
#         if !(is_linked_variable(var))
#             constraint_refs = JuMP.all_constraints(node_model, func, set)
#             for constraint_ref in constraint_refs
#                 constraint = JuMP.constraint_object(constraint_ref)
#                 new_constraint = _copy_constraint(constraint,reference_map)
#                 new_ref= JuMP.add_constraint(aggregate_model,new_constraint)
#                 agg_node.constraintmap[new_ref] = constraint_ref
#                 reference_map[constraint_ref] = new_ref
#             end
#         end
#     end
# end
