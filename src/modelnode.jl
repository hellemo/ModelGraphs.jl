##############################################################################
# Model Nodes
##############################################################################
"""
The ModelNode type

ModelNode()

Creates an empty ModelNode.  Does not add it to a graph.
"""
mutable struct ModelNode <: JuMP.AbstractModel
    #Meta data and index information for the hypergraph
    hypernode::HyperNode

    #The model
    model::JuMP.AbstractModel
    linkvariablemap::Dict{JuMP.AbstractVariableRef,AbstractLinkVariableRef}  #link variables this node uses

    #Solution Data
    variable_values::Dict{JuMP.AbstractVariableRef,Float64}               
    constraint_dual_values::Dict{JuMP.ConstraintRef,Float64}
    nl_constraint_dual_values::Dict{JuMP.NonlinearConstraintIndex,Float64}
end
#Constructor
ModelNode(hypernode::HyperNode) = ModelNode(hypernode,JuMP.Model(),Dict{JuMP.AbstractVariableRef,AbstractLinkVariableRef}(),Dict{MOI.VariableIndex,Float64}(),Dict{MOI.ConstraintIndex,Float64}(),Dict{JuMP.NonlinearConstraintIndex,Float64}())

JuMP.object_dictionary(m::ModelNode) = m.model.obj_dict
JuMP.variable_type(::ModelNode) = JuMP.VariableRef

function NestedHyperGraphs.add_node!(graph::AbstractModelGraph)
    hypergraph = gethypergraph(graph)
    hypernode = add_node!(hypergraph)
    node = ModelNode(hypernode)
    return node
end

function NestedHyperGraphs.add_node!(graph::AbstractModelGraph,m::AbstractModel)
    node = NestedHyperGraphs.add_node!(graph)
    set_model(node,m)
    return node
end

#Get node for a JuMP model if the model is set to a node
getnode(m::JuMP.Model) = m.ext[:modelnode]

#Get the corresponding node for a JuMP variable reference
function getnode(var::JuMP.VariableRef)
    if haskey(var.model.ext,:modelnode)
        return getnode(var.model)
    else
        error("variable $var does not belong to a modelnode.  If you're trying to create a linkconstraint, make sure
        the owning model has been set to a node.")
    end
end

function getnode(con::JuMP.ConstraintRef)
    if haskey(con.model.ext,:modelnode)
        return getnode(con.model)
    else
        error("constraint $con does not belong to a node")
    end
end

#Model Management
"Get the underlying JuMP model for a node"
getmodel(node::ModelNode) = node.model

"Get an underlying model variable"
#look for variables, if no variable, look for an attribute
function Base.getindex(node::ModelNode,symbol::Symbol)
    if haskey(node.model.obj_dict,symbol)
        return getmodel(node)[symbol]
    else
        return getattribute(node,symbol)
    end
end

function Base.setindex!(node::ModelNode,value::Any,symbol::Symbol)
    setattribute(node,symbol,value)
end


function JuMP.add_variable(node::ModelNode,  v::JuMP.AbstractVariable, name::String="")
    vref = JuMP.add_variable(getmodel(node),v,name)
    return vref
end

function JuMP.add_constraint(node::ModelNode,  con::JuMP.AbstractConstraint, name::String="")
    cref = JuMP.add_constraint(getmodel(node),con,name)          #also add to master model
    return cref
end

"""
JuMP.objective_function(node::ModelNode)

Get a node objective function.
"""
JuMP.objective_function(node::ModelNode) = JuMP.objective_function(getmodel(node))

"Get node objective value"
JuMP.objective_value(node::ModelNode) = JuMP.objective_value(node.model)

"""
linkconstraints(node::ModelNode)

Return a Dictionary of LinkConstraints for each graph the node is a member of
"""
# #TODO get incident edges to node and return those, or cache references on the node?
# function linkconstraints(node::ModelNode)
#     #links = Dict()
#     hypernode = gethypernode(node)
#
#     hyperedges = get_incident_edges(hypernode)  #This will return the hyper edges for each graph the hypernode is a part of
#
#     linkedges = getlinkedge.(hyperedges)
#     linkconstraints = getlinkconstraints.(linkedges)
#     # for (graph,refs) in node.linkconstraints
#     #     links[graph] = Vector{LinkConstraint}()
#     #     for ref in refs
#     #         push!(links[graph],JuMP.constraint_object(ref))
#     #     end
#     # end
#     return linkconstraints
# end

"""
linkconstraints(graph::AbstractModelGraph,node::ModelNode)

Return Array of LinkConstraints that cover the node
"""
# TODO
# function linkconstraints(graph::AbstractModelGraph,node::ModelNode)
#     links = []
#     for ref in node.linkconrefs[graph]
#         push!(links,LinkConstraint(ref))
#     end
#     return links
# end

########################################
# Get model node from other objects
########################################
"""
is_node_variable(node::ModelNode,var::AbstractJuMPScalar)

Check whether a JuMP variable belongs to a ModelNode
"""
is_node_variable(node::ModelNode,var::JuMP.AbstractVariableRef) = getmodel(node) == var.m   #checks whether a variable belongs to a node or edge

is_set_to_node(m::AbstractModel) = haskey(m.ext,:modelnode)                      #checks whether a model is assigned to a node

JuMP.num_variables(node::ModelNode) = JuMP.num_variables(getmodel(node))

########################################
#Get model nodes corresponding to models or variables
########################################
"""
NestedHyperGraphs.getnode(model::AbstractModel)

Get the ModelNode corresponding to a JuMP Model
"""
NestedHyperGraphs.getnode(m::AbstractModel) = is_set_to_node(m) ? m.ext[:modelnode] : throw(error("Only node models have associated graph nodes"))

"""
StructureGraphs.getnode(model::AbstractModel)

Get the ModelNode corresponding to a JuMP Variable
"""
NestedHyperGraphs.getnode(var::JuMP.AbstractVariableRef) = JuMP.owner_model(var).ext[:node]

"""
set_model(node::ModelNode,m::AbstractModel)

Set the model on a node.  This will delete any link-constraints the node is currently part of
"""
function set_model(node::ModelNode,m::JuMP.AbstractModel;preserve_links = false)
    !(is_set_to_node(m) && getmodel(node) == m) || error("Model $m is already asigned to another node")
    # TODO
    # BREAK LINKS FOR NOW
    # Check for the same variable containers to attach link constraints to
    # If it already had a model, delete all the link constraints corresponding to that model
    # if hasmodel(node)
    #     for (graph,constraints) in getlinkconstraints(node)
    #         local_link_cons = constraints
    #         graph_links = getlinkconstraints(graph)
    #         filter!(c -> !(c in local_link_cons), graph_links)  #filter out local link constraints
    #         node.link_data = NodeLinkData()   #reset the local node or edge link data
    #     end
    # end
    node.model = m
    m.ext[:modelode] = node
end

#TODO
#set a model with the same variable names and dimensions as the old model on the node.
#This will not break link constraints by default but will make sure they match the old model
#switch out variables in any connected linkconstraints
#throw warnings if link constraints break
function resetmodel(node::ModelNode,m::JuMP.AbstractModel)
    #reassign the model
    node.model = m
end

function string(node::ModelNode)
    "Model Node w/ $(JuMP.num_variables(node)) Variable(s)"
end
print(io::IO,node::ModelNode) = print(io, string(node))
show(io::IO,node::ModelNode) = print(io,node)

# TODO
# clearmodel(node::ModelNode) = nodeoredge.attributes[:model] = nothing  #need to update link constraints

getnodevariable(node::ModelNode,index::Integer) = JuMP.VariableRef(getmodel(node),index)

JuMP.all_variables(node::ModelNode) = JuMP.all_variables(getmodel(node))
nodevalue(var::JuMP.VariableRef) = getnode(var).variable_values[var.index]  #TODO #Get values of JuMP expressions
nodedual(con_ref::JuMP.ConstraintRef{JuMP.Model,MOI.ConstraintIndex}) = getnode(con).constraint_dual_values[con.index]
nodedual(con_ref::JuMP.ConstraintRef{JuMP.Model,JuMP.NonlinearConstraintIndex}) = getnode(con).nl_constraint_dual_values[con.index]
