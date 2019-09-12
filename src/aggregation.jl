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
getlinkconstraints(m::JuMP.Model) = is_aggregate_model(m) && getaggregationinfo(m).linkconstraints
getlinkvariables(m::JuMP.Model) = is_aggregate_model(m) && getaggregationinfo(m).linkvariables
getNLlinkconstraints(m::JuMP.Model) = is_aggregate_model(m) && getaggregationinfo(m).NLlinkconstraints
NHG.getnodes(m::JuMP.Model) = is_aggregate_model(m) && getaggregationinfo(m).nodes

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
function aggregate(modelgraph::ModelGraph)
    aggregate_model = AggregateModel()
    reference_map = AggregationMap(aggregate_model)

    master_reference_map = _add_to_aggregate_model!(aggregate_model,getmastermodel(modelgraph),reference_map)
    #merge!(reference_map,master_reference_map)

    #COPY NODE MODELS INTO AGGREGATED MODEL
    has_nonlinear_objective = false                     #check if any nodes have nonlinear objectives
    for modelnode in getnodes(modelgraph)               #for each node in the model graph
        node_model = getmodel(modelnode)
        #Need to pass master reference so we use those variables instead of creating new ones
        _add_to_aggregate_model!(aggregate_model,node_model,reference_map)  #updates jump_graph_model,the jump_node, and the ref_map
        #merge!(reference_map,node_reference_map)   #Update the reference_map

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
        error("NL graph objective not yet supported")
        # dgraph = JuMP.NLPEvaluator(modelgraph)
        # MOI.initialize(dgraph,[:ExprGraph])
        # graph_obj = MOI.objective_expr(dgraph)
        # _splice_nonlinear_variables!(graph_obj,reference_map)  #_splice_nonlinear_variables!(node_obj,var_maps[node])
        # JuMP.set_NL_objective(aggregate_model,JuMP.objective_sense(modelgraph,graph_obj))
    else
        _set_node_objectives!(modelgraph,aggregate_model,reference_map,has_nonlinear_objective)          #set to the sum of the node objectives
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

#Aggregate a graph based on a model partition.  Return a new ModelGraph with possible subgraphs (If it was passed a recursive partition)
function aggregate(graph::ModelGraph,hyperpartition::HyperPartition)
    println("Building Aggregate Model Graph using HyperPartition")

    #Create New ModelGraphs
    parent_dict = Dict()
    for parent in hyperpartition.parents
        new_model_graph = ModelGraph()
        parent_dict[parent] = new_model_graph
    end

    top_model_graph = parent_dict[hyperpartition.parents[1]]
    reference_map = AggregationMap(top_model_graph)  #old model graph => new modelgraph


    #BOTTOM LEVEL NODES
    #Aggregate subgraphs to create bottom level nodes
    submodelgraphs = []
    for partition in hyperpartition.partitions
        hypergraph = partition.hypergraph
        submodelgraph = ModelGraphs.create_sub_modelgraph(graph,hypergraph)
        push!(submodelgraphs,submodelgraph)

        aggregate_model,agg_ref_map = aggregate(submodelgraph)

        merge!(reference_map,agg_ref_map)

        parent_graph = parent_dict[partition.parent]
        aggregate_node = add_node!(parent_graph)
        set_model(aggregate_node,aggregate_model)
    end


    # #Now add shared nodes and shared edges
    for parent in hyperpartition.parents
        shared_nodes = parent.sharednodes     #Could be linkconstraints, shared variables, shared models, or pairs
        shared_edges = parent.sharededges

        parent_mg = parent_dict[parent]

        #LINK VARIABLES
        # master = aggregate(shared_nodes) #get linkvariables from shared nodes
        # set_master(parent_mg,master)
        master = Model()
        for shared_node in shared_nodes
            error("Shared nodes not supported yet")
            #identify edges here and figure out which link variables to make
        end
        parent_mg.mastermodel = master

        #LINK CONSTRAINTS
        for shared_edge in shared_edges
            linkedge = findlinkedge(graph,shared_edge)
            for linkconstraintref in linkedge.linkconstraints
                linkconstraint = LinkConstraint(linkconstraintref)
                new_con = _copy_constraint(linkconstraint,reference_map)
                JuMP.add_constraint(parent_mg,new_con)  #this is a link constraint
            end
        end

        if parent.parent != nothing
            parent_subgraph = parent_dict[parent.parent]
            add_subgraph!(subgraph,new_model_graph)
        end
    end
    return top_model_graph,reference_map
end

#Aggregate the subgraphs of a modelgrap where n_levels corresponds to how many levels remain, 0 means no subgraphs
function aggregate(graph::ModelGraph,n_levels::Int64)
    new_model_graph = ModelGraph()
    return new_model_graph
end


#previously _buildnodemodel!
function _add_to_aggregate_model!(aggregate_model::JuMP.Model,node_model::JuMP.Model,aggregation_map::AggregationMap)          #jump_node

    agg_node = add_aggregated_node!(aggregate_model)

    if JuMP.mode(node_model) == JuMP.DIRECT
        error("Cannot copy a node model in `DIRECT` mode. Use the `Model` ",
              "constructor instead of the `direct_model` constructor to be ",
              "able to aggregate into a new JuMP Model.")
    end

    #reference_map = GraphReferenceMap(m,MOIU.IndexMap())
    reference_map = AggregationMap(aggregate_model)


    constraint_types = JuMP.list_of_constraint_types(node_model)
    #COPY VARIABLES
    for var in JuMP.all_variables(node_model)
        if is_linked_variable(var)                                       #if the variable is actually a link variable, we don't need to make a new one
            reference_map[var] = aggregation_map[getlinkvariable(var)]                    #get the master variable
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
            _splice_nonlinear_variables!(node_obj,node_model,reference_map)  #_splice_nonlinear_variables!(node_obj,var_maps[node])
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
    new_obj = Expr(:call,:*,:($sense),new_obj)
    return new_obj
end

function set_sum_of_objectives(graph::ModelGraph)
end
