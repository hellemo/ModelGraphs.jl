##############################################################################
# Model Nodes
##############################################################################
"""
The ModelNode type

ModelNode()

Creates an empty ModelNode.  Does not add it to a graph.
"""
mutable struct ModelNode <: AbstractModelNode
    structurenode::StructureGraphs.StructureNode
    model::JuMP.AbstractModel

    # NOTE: Thinking whether we store linkconstraint references.  Depends how often this would have to be accessed.
    #linkconstraints::Dict{AbstractModelGraph,Vector{ConstraintRef}}
    #linkconstraints::DefaultDict{AbstractModelGraph, Vector{GraphConstraintRef}}(GraphConstraintRef[])
end
#Constructor
ModelNode() = ModelNode(StructureGraphs.StructureNode(),JuMP.Model())#,Dict{AbstractModelGraph,Vector{ConstraintRef}}())
StructureGraphs.create_node(graph::AbstractModelGraph) = ModelNode()
StructureGraphs.getstructurenode(node::ModelNode) = node.structurenode

"""
addnode!(graph::AbstractModelGraph)

Add a ModelNode to a ModelGraph.
"""
function StructureGraphs.add_node!(graph::AbstractModelGraph,m::AbstractModel)
    node = StructureGraphs.add_node!(graph)
    setmodel(node,m)
    return node
end

#Get node for a JuMP model if the model is set to a node
StructureGraphs.getnode(m::JuMP.Model) = m.ext[:node]

#Get the corresponding node for a JuMP variable reference
function StructureGraphs.getnode(var::JuMP.VariableRef)
    if haskey(var.model.ext,:node)
        return getnode(var.model)
    else
        error("variable $var does not belong to a node.  If you're trying to create a linkconstraint, make sure
        the owning model has been set to a node.")
    end
end

#Model Management
"Get the underlying JuMP model for a node"
getmodel(node::ModelNode) = node.model

"Get an underlying model variable"
#look for variables, if no variable, look for an attribute
function getindex(node::ModelNode,symbol::Symbol)
    if haskey(node.model.obj_dict,symbol)
        return getmodel(node)[symbol]
    else
        return getattribute(node,symbol)
    end

end
#getindex(node::ModelNode,sym::Symbol) = getmodel(node)[sym]  #get a variable on a node

"""
getobjective(node::ModelNode)

Get a node objective.
"""
#TODO Update to JuMP 0.19
# JuMP.getobjective(node::ModelNode) = getobjective(node.model)

"Get node objective value"
JuMP.getobjectivevalue(node::ModelNode) = getobjectivevalue(node.model)

"""
getlinkconstraints(node::ModelNode)

Return a Dictionary of LinkConstraints for each graph the node is a member of
"""
#TODO get incident edges to node and return those, or cache references on the node?
# function getlinkconstraints(node::ModelNode)
#     links = Dict()
#     for (graph,refs) in node.linkconstraints
#         links[graph] = Vector{LinkConstraint}()
#         for ref in refs
#             push!(links[graph],JuMP.constraint_object(ref))
#         end
#     end
#     return links
# end

"""
getlinkconstraints(graph::AbstractModelGraph,node::ModelNode)

Return Array of LinkConstraints that cover the node
"""
# TODO
# function getlinkconstraints(graph::AbstractModelGraph,node::ModelNode)
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
is_set_to_node(m::AbstractModel) = haskey(m.ext,:node)                      #checks whether a model is assigned to a node
JuMP.num_variables(node::ModelNode) = JuMP.num_variables(getmodel(node))

########################################
#Get model nodes corresponding to models or variables
########################################
"""
getnode(model::AbstractModel)

Get the ModelNode corresponding to a JuMP Model
"""
StructureGraphs.getnode(m::AbstractModel) = is_set_to_node(m) ? m.ext[:node] : throw(error("Only node models have associated graph nodes"))

"""
getnode(model::AbstractModel)

Get the ModelNode corresponding to a JuMP Variable
"""
StructureGraphs.getnode(var::JuMP.AbstractVariableRef) = JuMP.owner_model(var).ext[:node]

"""
setmodel(node::ModelNode,m::AbstractModel)

Set the model on a node.  This will delete any link-constraints the node is currently part of
"""
function setmodel(node::ModelNode,m::AbstractModel;preserve_links = false)
    !(is_set_to_node(m) && getmodel(node) == m) || error("the model is already asigned to another node")
    #TODO
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
    m.ext[:node] = node
end

#TODO
#set a model with the same variable names and dimensions as the old model on the node.
#This will not break link constraints by default but will make sure they match the old model
#switch out variables in any connected linkconstraints
#throw warnings if link constraints break
function resetmodel(node::ModelNode,m::AbstractModel)
    #reassign the model
    node.model = m
end

function string(node::ModelNode)
    "Model Node w/ $(JuMP.num_variables(node)) Variable(s)"
end
print(io::IO,node::ModelNode) = print(io, string(node))
show(io::IO,node::ModelNode) = print(io,node)

# TODO
# removemodel(node::ModelNode) = nodeoredge.attributes[:model] = nothing  #need to update link constraints

# TODO Rewrite
# getnodevariable(node::ModelNode,index::Integer) = Variable(getmodel(node),index)
# getnodevariable(node::ModelNode)

#TODO Rewrite for new JuMP v0.19
# This effectively return a mapping of symbols to different JuMP containers. This can be done better.
# function getnodevariablemap(node::ModelNode)
#     node_map = Dict()
#     node_model = getmodel(node)
#     for key in keys(node_model.objDict)  #this contains both variable and constraint references
#         if isa(node_model.objDict[key],Union{JuMP.JuMPArray{AbstractJuMPScalar},Array{AbstractJuMPScalar}})     #if the JuMP variable is an array or a JuMPArray
#             vars = node_model.objDict[key]
#             node_map[key] = vars
#         #reproduce the same mapping in a dictionary
#         elseif isa(node_model.objDict[key],JuMP.JuMPDict)
#             tdict = node_model.objDict[key].tupledict  #get the tupledict
#             d_tmp = Dict()
#             for dkey in keys(tdict)
#                 d_tmp[dkey] = var_map[linearindex(tdict[dkey])]
#             end
#             node_map[key] = d_tmp
#
#         elseif isa(node_model.objDict[key],JuMP.AbstractJuMPScalar) #else it's a single variable
#             node_map[key] = node_model.objDict[key]
#         # else #objDict also has contraints!
#         #     error("Did not recognize the type of a JuMP variable $(node_model.objDict[key])")
#         end
#     end
#     return node_map
# end

# getlinkreferences(node::ModelNode) = node.linkconrefs
# getlinkreferences(graph::AbstractModelGraph,node::ModelNode) = node.linkconrefs[graph]
