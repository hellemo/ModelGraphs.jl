"""
    AggregateGraph
    Extension Data attached to an Aggregated ModelGraph.  The AggregateGraph retains a reference to the original ModelGraph topology with references to
    graphvariables, graphconstraints, and linkconstraints.
"""
#Aggregation Data
#After aggregating, we throw away the topology, but we keep references to what was a linkvariable and linkconstraint
mutable struct AggregationInfo # <: AbstractModelGraph
    #hypergraph::NestedHyperGraph
    nodes::Vector{NodeReference}
    linkvariables::Vector{VariableRef}
    linkconstraints::Vector{ConstraintRef}
    NLlinkconstraint::Vector{ConstraintRef}
end
AggregateInfo() = AggregateModel()
gethypergraph(graph::AggregateGraph) = graph.hypergraph

mutable struct NodeReference #<: AbstractModelNode
    hypernode::HyperNode
    obj_dict::Dict{Symbol,Any}
    variablemap::Dict{JuMP.VariableRef,JuMP.VariableRef}
    constraintmap::Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}
    nl_constraintmap::Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}
    objective::Union{JuMP.AbstractJuMPScalar,Expr}
end

zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef}))

#Construct a structured model, but roll it all into one JuMP model (this is how we solve with JuMP accessible solvers)
function AggregateModel()
    m = JuMP.Model()
    m.ext[:AggregationInfo] = AggregationInfo()
    return m
end

is_aggregate_model(m::JuMP.Model) = haskey(m.ext,:AggregationInfo) ? true : false  #check if the model is a graph model
assert_aggregate_model(m::JuMP.Model) = @assert is_aggregate_model(m)
#Define all of the JuMP model extension functions
get_aggregation_info(m::JuMP.Model) = haskey(m.ext, :AggregationInfo) ? m.ext[:AggregationInfo] : error("Model is not a graph model")


JuMP.objective_function(node::NodeReference) = node.objective

getnodevariables(node::NodeReference) = node.variablelist
getnodevariable(node::NodeReference,index::Integer) = node.variablelist[index]
getnodeconstraints(node::NodeReference) = node.constraintlist
JuMP.num_variables(node::NodeReference) = length(node.variablelist)
#get node variables using a symbol lookup
getindex(node::NodeReference,s::Symbol) = node.obj_dict[s]

#get all of the link constraints from a JuMP model
#Retrieves constraint indices that are link constraints
function get_link_constraints(m::JuMP.Model)
    assert_aggregate_model(m)
    agg_info = get_aggregation_info(m)
    return agg_info.linkconstraints
end

"""
    AggregationMap
    Mapping between variable and constraint reference of a node model to the Aggregated Model.
    The reference of the aggregated model can be obtained by indexing the map with the reference of the corresponding reference of the original node model.
"""
struct AggregationMap
    aggregate_model::JuMP.Model   #An aggregate model
    varmap::Dict{JuMP.VariableRef,JuMP.VariableRef}
    conmap::Dict{JuMP.ConstraintRef,JuMP.ConstraintRef}
end
function Base.getindex(reference_map::AggregationMap, vref::JuMP.VariableRef)  #reference_map[node_var] --> aggregated_copy_var
    #return JuMP.VariableRef(reference_map.model,reference_map.index_map[JuMP.index(vref)])
    return reference_map.varmap[vref]
end
function Base.getindex(reference_map::AggregationMap, cref::JuMP.ConstraintRef)
    #return JuMP.ConstraintRef(reference_map.model,reference_map.index_map[JuMP.index(cref)],cref.shape)
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


function aggregate(graph::ModelGraph;levels = 0)
    node = Modelnode()

    return node
end

function aggregate(graph::ModelGraph,partition::ModelPartition)
end

function aggregate(graph::ModelGraph,nodes::Vector{ModelNode})
end

#Create an aggregated ModelGraph that can be used by decomposition algorithms
function create_aggregate_graph(model_graph::ModelGraph,partition_data::PartitionData)
    println("Building Aggregated Model Graph")
    new_model_graph = ModelGraph()  #The new aggregated ModelGraph

    partitions = partition_data.partitions
    shared_entities = partition_data.shared_entities #Could be linkconstraints, shared variables, shared models, or pairs

    #What kind of reference map do we need?  multiple?
    variable_map = Dict{JuMP.VariableRef,JuMP.VariableRef}()
    n_partitions = length(partitions)
    #Aggregate Partitions
    for i = 1:n_partitions  #create aggregate model for each partition
        part = partitions[i]

        #Get LinkingEdges from subgraphs if there are any.
        local_shared_entities = partition_data.partition_entities[i]

        aggregate_model,agg_ref_map = create_aggregate_model(model_graph,part,local_shared_entities)

        #Update VariableMap
        merge!(variable_map,agg_ref_map.varmap)

        aggregate_node = add_node!(new_model_graph)
        setmodel(aggregate_node,aggregate_model)

    end

    #NEW LINK CONSTRAINTS
    for entity in shared_entities #Could be e.g. LinkConstraints
        add_shared_entity!(new_model_graph,entity,variable_map)
    end

    return new_model_graph

end

function add_shared_entity!(graph::ModelGraph,edge::LinkingEdge,variable_map::Dict)
    link_constraints = [LinkConstraint(link_ref) for link_ref in edge.linkconstraints]
    for link_constraint in link_constraints
        copy_constraint = AlgebraicGraphs._copy_constraint(link_constraint,variable_map)
        JuMP.add_constraint(graph.linkmodel,copy_constraint)
    end
end

function add_shared_entity!(graph::ModelGraph,graphvariable::JuMP.AbstractVariableRef,ref_map::GraphReferenceMap)
end

#Build up an aggregate model given a set of nodes.
#function create_aggregate_model(model_graph::ModelGraph,nodes::Vector{ModelNode},link_edges::Vector{LinkingEdge})
function create_aggregate_model(model_graph::ModelGraph,nodes::Vector,link_edges::Vector)
    #local_links, cross_links = _get_local_and_cross_links(model_graph,nodes)

    #Get corresponding ModelNodes and LinkConstraints for given indices
    #Extract LinkConstraints from the Edges
    link_constraints = []
    for edge in link_edges
        for link_ref in edge.linkconstraints
            linkcon = LinkConstraint(link_ref)
            push!(link_constraints,linkcon)
        end
    end

    aggregate_model =  AlgebraicGraphs.JuMPGraphModel()         #Use a JuMPGraphModel so we can track the internal structure
    jump_graph = getgraph(aggregate_model)

    reference_map = GraphReferenceMap(aggregate_model,Dict(),Dict())
    #variable_map = Dict()

    has_nonlinear_objective = false

    #COPY NODE MODELS INTO AGGREGATE MODEL
    for model_node in nodes  #for each node in the model graph
        nodeindex = getindex(model_graph,model_node)
        jump_node = add_node!(aggregate_model,index = nodeindex)  #add at the same index

        node_reference_map = _buildnodemodel!(aggregate_model,jump_node,model_node)

        merge!(reference_map,node_reference_map)

        node_model = getmodel(model_node)

        if has_nonlinear_objective != true
            has_nonlinear_objective = _has_nonlinear_obj(node_model)
        end

        #TODO Get nonlinear object data to work
        #COPY OBJECT DATA (JUMP CONTAINERS).  I don't really need this for this.  It would be nice for Aggregation though.
        for (name, value) in JuMP.object_dictionary(node_model)
            #jump_node.obj_dict[name] = reference_map[value]
            jump_node.obj_dict[name] = getindex.(reference_map, value)
        end
    end

    #LOCAL LINK CONSTRAINTS
    #TODO: Typing is causing issues here
    for linkconstraint in link_constraints
        new_constraint = _copy_constraint(linkconstraint,reference_map)
        JuMP.add_constraint(aggregate_model,new_constraint)
    end

    #NODE OBJECTIVE
    if !(has_nonlinear_objective)
        graph_obj = sum(JuMP.objective_function(jump_node) for jump_node in getnodes(jump_graph)) #NOTE: All of the node objectives were converted to Minimize (MOI.OptimizationSense(0))
        JuMP.set_objective(aggregate_model,MOI.OptimizationSense(0),graph_obj)
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
            _splice_nonlinear_variables!(node_obj,ref_map)  #_splice_nonlinear_variables!(node_obj,var_maps[node])

            node_obj = Expr(:call,:*,:($sense),node_obj)
            graph_obj = Expr(:call,:+,graph_obj,node_obj)
        end
        JuMP.set_NL_objective(aggregate_model, MOI.OptimizationSense(0), graph_obj)
    end

    #Need to return multiple reference maps

    return aggregate_model,reference_map

end

#Build up an aggregate model given a set of edges
function create_aggregate_model(model_graph::ModelGraph,edges::Vector{LinkingEdge},shared_variables::Vector{GraphVariableRef})
end

#Function to build a node model for a flat graph model
function _add_to_aggregate_model!(m::JuMP.Model,model_node::ModelNode)  #jump_node

    #m = aggregate_model
    node_reference = add_ref_node!(m)

    node_model = getmodel(model_node)

    if JuMP.mode(node_model) == JuMP.DIRECT
        error("Cannot copy a node model in `DIRECT` mode. Use the `Model` ",
              "constructor instead of the `direct_model` constructor to be ",
              "able to aggregate into a new JuMP Model.")
    end

    #reference_map = GraphReferenceMap(m,MOIU.IndexMap())
    reference_map = AggregationMap(m)

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

        node_reference.variablemap[new_x] = var
        #jump_node.variablemap[new_x] = var
        #push!(jump_node.variablelist,new_x)
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


            new_ref= JuMP.add_constraint(m,new_constraint)
            #push!(jump_node.constraintlist,ref)
            #jump_node.constraintmap[new_ref] = constraint_ref
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
            new_nl_constraint = JuMP.add_NL_constraint(m,expr)      #raw expression input for non-linear constraint
            constraint_ref = JuMP.ConstraintRef(node_model,JuMP.NonlinearConstraintIndex(i),new_nl_constraint.shape)
            #push!(jump_node.nl_constraints,constraint_ref)
            #jump_node.nl_constraintmap[new_nl_constraint] = constraint_ref
            node_reference.nl_constraintmap[new_nl_constraint] = constraint_ref
            reference_map[constraint_ref] = new_nl_constraint
        end
    end

    #OBJECTIVE FUNCTION
    if !(_has_nonlinear_obj(node_model))
        #AFFINE OR QUADTRATIC OBJECTIVE
        new_objective = _copy_objective(node_model,reference_map)
        #jump_node.objective = new_objective
        node_reference.objective = new_objective
    else
        #NONLINEAR OBJECTIVE
        if !nlp_initialized
            d = JuMP.NLPEvaluator(node_model)           #Get the NLP evaluator object.  Initialize the expression graph
            MOI.initialize(d,[:ExprGraph])
        end
        new_obj = _copy_nl_objective(d,variablemap)
        #jump_node.objective = new_obj
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
