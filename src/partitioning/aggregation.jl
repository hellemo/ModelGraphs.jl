#Create an aggregated ModelGraph that can be used by decomposition algorithms
function create_aggregate_graph(model_graph::ModelGraph,partition_data::PartitionData)
    new_model_graph = ModelGraph()  #The new aggregated ModelGraph

    partitions = partition_data.partitions
    shared_entities = partition_data.shared_entities #Could be linkconstraints, shared variables, shared models, or pairs

    #What kind of reference map do we need?  multiple?
    #reference_map = AlgebraicGraphs.GraphReferenceMap()  What should this have?
    variable_map = Dict()

    n_partitions = length(partitions)
    #Aggregate Partitions
    for i = 1:n_partitions  #create aggregate model for each partition
        part = partitions[i]

        #Need to get LinkingEdges from subgraphs
        local_shared_entities = partition_data.partition_entities[i]

        aggregate_model,agg_ref_map = create_aggregate_model(model_graph,part,local_shared_entities)

        #Update ReferenceMap
        #merge!(reference_map,agg_ref_map)
        #TODO: This is not correct.  Values get overridden
        merge!(variable_map,agg_ref_map.varmap)

        aggregate_node = add_node!(new_model_graph)
        setmodel(aggregate_node,aggregate_model)

    end

    #NEW LINK CONSTRAINTS
    for entity in shared_entities #Could be e.g. LinkConstraints
        add_shared_entity!(new_model_graph,entity,variable_map)  #e.g. copy the LinkConstraints onto a new ModelGraph
    end

    return new_model_graph

end

function add_shared_entity!(graph::ModelGraph,edge::LinkingEdge,variable_map::Dict)
    link_constraints = [LinkConstraint(link_ref) for link_ref in edge.linkconstraints]
    for link_constraint in link_constraints
        copy_constraint = AlgebraicGraphs._copy_constraint(link_constraint,variable_map)
        JuMP.add_constraint(graph.link_model,copy_constraint)
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

    reference_map = AlgebraicGraphs.GraphReferenceMap(aggregate_model,Dict(),Dict())
    #variable_map = Dict()

    has_nonlinear_objective = false

    #COPY NODE MODELS INTO AGGREGATE MODEL
    for model_node in nodes  #for each node in the model graph
        nodeindex = getindex(model_graph,model_node)
        jump_node = add_node!(aggregate_model,index = nodeindex)  #add at the same index

        node_reference_map = AlgebraicGraphs._buildnodemodel!(aggregate_model,jump_node,model_node)


        merge!(reference_map,node_reference_map)
        # m,var_map = _buildnodemodel!(aggregate_model,jump_node,model_node)  #build model nodes into new jump node.  var_map is a mapping of model node variable indices to new jump node variables
        # var_maps[model_node] = var_map

        node_model = getmodel(model_node)

        if has_nonlinear_objective != true
            has_nonlinear_objective = AlgebraicGraphs._has_nonlinear_obj(node_model)
        end
    end

    #LOCAL LINK CONSTRAINTS
    #TODO: Typing is causing issues here
    #BUG: Some nodes in the link constraint are NOT in the set of nodes.  Identifying shared entities is either wrong, or the mapping is incorrect.
    for linkconstraint in link_constraints
        new_constraint = AlgebraicGraphs._copy_constraint(linkconstraint,reference_map)
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


#Given a model graph and a set of node partitions, create a new aggregated model graph
# #NOTE: I'll have an option to transform linkconstraints into shared variables
# function create_lifted_model_graph(model_graph::ModelGraph,partitions::Vector{Vector{ModelNode}})
#     new_model_graph = ModelGraph()
#
#     aggregate_models = []
#     all_cross_links = []
#     variable_mapping = Dict()
#     all_var_maps = Dict()
#     for partition in partitions  #create aggregate model for each partition
#         aggregate_model,cross_links,var_maps = create_aggregate_model(model_graph,partition)
#         push!(aggregate_models,aggregate_model)
#         append!(all_cross_links,cross_links)
#         merge!(all_var_maps,var_maps)
#     end
#     all_cross_links = unique(all_cross_links)  #remove duplicate cross links
#
#     for agg_model in aggregate_models
#         aggregate_node = add_node!(new_model_graph)
#         setmodel(aggregate_node,agg_model)
#     end
#
#     sub_nodes = collectnodes(new_model_graph)
#     #GLOBAL LINK CONSTRAINTS.  Re-add link constraints to aggregated model nodes
#     master_node = add_node!(new_model_graph,Model())
#
#     for linkconstraint in all_cross_links
#         linear_terms = []
#         for terms in linearterms(linkconstraint.terms)
#             push!(linear_terms,terms)
#         end
#
#         #Get references to variables in the aggregated models
#         t_new = []
#         for i = 1:length(linear_terms)
#             coeff = linear_terms[i][1]
#             var = linear_terms[i][2]
#             model_node = getnode(var)                #the original model node
#             var_index = JuMP.linearindex(var)        #variable index in the model node
#             var_map = all_var_maps[model_node]       #model node variable map {index => aggregate_variable}
#             agg_var = var_map[var_index]
#             push!(t_new,(coeff,agg_var))
#         end
#
#         #For simple links, create a single lifted variable and 2 linkconstraints
#         if length(t_new) == 2 && (linkconstraint.lb == linkconstraint.ub) && (linear_terms[1][1]) == -linear_terms[2][1]
#             lifted_var = @variable(getmodel(master_node))
#             for term in t_new
#                 var = term[2]
#                 @linkconstraint(new_model_graph,lifted_var == var)
#             end
# #        end
#
#         else  #else create multiple duplicates and lift the constraints into the master
#             lift_terms = []
#             for term in t_new
#                 coeff = term[1]
#                 var = term[2]
#                 lifted_var = @variable(getmodel(master_node) , start = JuMP.getvalue(var))
#                 @linkconstraint(new_model_graph,lifted_var == var)
#                 push!(lift_terms,(coeff,lifted_var))
#             end
#             @constraint(getmodel(master_node), linkconstraint.lb <= sum(lift_terms[i][1]*lift_terms[i][2] for i = 1:length(lift_terms)) + linkconstraint.terms.constant <= linkconstraint.ub)
#         end
#
#         #For simple links, create a single lifted variable and 2 linkconstraints
#         @linkconstraint(new_model_graph, linkconstraint.lb <= sum(t_new[i][1]*t_new[i][2] for i = 1:length(t_new)) + linkconstraint.terms.constant <= linkconstraint.ub)
#     end
#     return new_model_graph , master_node, sub_nodes
# end
