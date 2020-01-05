##############################################################################
# ModelGraph
##############################################################################
"""
ModelGraph()

The ModelGraph Type.  Represents a graph containing models (nodes) and the linkconstraints (edges) between them.
"""
mutable struct ModelGraph <: AbstractModelGraph

    masternode::ModelNode                       #Master node in a graph.  Can contain linkvariables

    #Topology
    modelnodes::Vector{ModelNode}                #Local model nodes
    linkedges::Vector{LinkEdge}                  #Local link edges.  These can also connect nodes across subgraphs
    node_idx_map::Dict{ModelNode,Int64}          #Local map of model nodes to indices
    edge_idx_map::Dict{LinkEdge,Int64}           #Local map of link edges indices
    subgraphs::Vector{AbstractModelGraph}        #Subgraphs contained in the model graph

    #graphindex::Int64
    linkedge_map::OrderedDict{Set,LinkEdge}      #Sets of vertices map to a linkedge

    #Link variables
    linkvariables::OrderedDict{Int64,AbstractLinkVariableRef}                                 #Link Variable references to master node variables
    child_linkvariable_map::Dict{AbstractLinkVariableRef,Vector{JuMP.AbstractVariableRef}}    #Map of link variables in master model to corresponding variables in child ModelNodes.
    parent_linkvariable_map::Dict{JuMP.AbstractVariableRef,AbstractLinkVariableRef}           #Map of graph link variables to parent graph link variables
    linkvariable_names::Dict{Int64,String}

    #Link constraints
    linkconstraints::OrderedDict{Int64,AbstractLinkConstraint}                     #Link constraint.  Defined over variables in ModelNodes.
    linkeqconstraints::OrderedDict{Int64,AbstractLinkConstraint}
    linkineqconstraints::OrderedDict{Int64,AbstractLinkConstraint}
    linkconstraint_names::Dict{Int64,String}

    #Objective
    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar

    #Optimizer
    optimizer::Union{JuMP.OptimizerFactory,AbstractGraphOptimizer,Nothing}

    #Object indices
    linkvariable_index::Int
    linkeqconstraint_index::Int           #keep track of constraint indices
    linkineqconstraint_index::Int
    linkconstraint_index::Int

    #Model Information
    obj_dict::Dict{Symbol,Any}

    #TODO Nonlinear Link Constraints using NLP Data
    nlp_data::Union{Nothing,JuMP._NLPData}

    #Constructor
    function ModelGraph()
        modelgraph = new(ModelNode(),
                    Vector{ModelNode}(),
                    Vector{LinkEdge}(),
                    Dict{ModelNode,Int64}(),
                    Dict{LinkEdge,Int64}(),
                    Vector{AbstractModelGraph}(),
                    OrderedDict{Set,LinkEdge}(),
                    OrderedDict{Int,AbstractLinkVariableRef}(),
                    Dict{AbstractLinkVariableRef,Vector{JuMP.AbstractVariableRef}}(),
                    Dict{JuMP.AbstractVariableRef,AbstractLinkVariableRef}(),
                    Dict{Int,String}(),
                    OrderedDict{Int, AbstractLinkConstraint}(),
                    OrderedDict{Int, AbstractLinkConstraint}(),
                    OrderedDict{Int, AbstractLinkConstraint}(),
                    OrderedDict{Int64,String}(),
                    MOI.FEASIBILITY_SENSE,
                    zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef}),
                    nothing,
                    0,
                    0,
                    0,
                    0,
                    Dict{Symbol,Any}(),
                    nothing
                    )

        modelgraph.masternode.ext[:modelgraph] = true
        modelgraph.node_idx_map[modelgraph.masternode] = 0

        return modelgraph
    end
end


########################################################
# ModelGraph Interface
########################################################
#################
#Subgraphs
#################
function add_subgraph!(graph::ModelGraph,subgraph::ModelGraph)
    push!(graph.subgraphs,subgraph)
    return graph
end
#Recursively grab subgraphs in the given modelgraph
getsubgraphs(modelgraph::ModelGraph) = modelgraph.subgraphs

function all_subgraphs(modelgraph::ModelGraph)
    subgraphs = modelgraph.subgraphs
    for subgraph in subgraphs
        subgraphs = [subgraphs;all_subgraphs(subgraph)]
        #append!(subgraphs,all_subgraphs(subgraph))
    end
    return subgraphs
end

#################
#ModelNodes
#################
function add_node!(graph::ModelGraph)
    modelnode = ModelNode()
    push!(graph.modelnodes,modelnode)
    graph.node_idx_map[modelnode] = length(graph.modelnodes)
    return modelnode
end

function add_node!(graph::ModelGraph,m::JuMP.Model)
    node = add_node!(graph)
    set_model(node,m)
    return node
end

function add_node!(graph::ModelGraph,modelnode::ModelNode)
    push!(graph.modelnodes,modelnode)
    graph.node_idx_map[modelnode] = length(graph.modelnodes)
    return modelnode
end

getnodes(graph::ModelGraph) = graph.modelnodes
getnode(graph::ModelGraph,index::Int64) = graph.modelnodes[index]

#Recursively collect nodes in a modelgraph from each of its subgraphs
function all_nodes(graph::ModelGraph)
    nodes = graph.modelnodes
    for subgraph in graph.subgraphs
        nodes = [nodes;all_nodes(subgraph)]
    end
    return nodes
end

#Find a node from recursive collection of modelgraph nodes.
function find_node(graph::ModelGraph,index::Int64)
    nodes = getnodes(graph)
    return nodes[index]
end

function Base.getindex(graph::ModelGraph,node::ModelNode)
    return graph.node_idx_map[node]
end

#################
#LinkEdges
#################
function add_link_edge!(graph::ModelGraph,modelnodes::Vector{ModelNode})
    #node_indices = Set([graph.node_idx_map[node] for node in modelnodes])
    #Check for existing linkedge
    key = Set(modelnodes)
    if haskey(graph.linkedge_map,key)
        linkedge = graph.linkedge_map[key]
    else
        linkedge = LinkEdge(modelnodes)
        push!(graph.linkedges,linkedge)
    end
    n_links = length(graph.linkedges)
    idx = n_links + 1
    graph.linkedge_map[linkedge.nodes] = linkedge
    graph.edge_idx_map[linkedge] = idx
    return linkedge
end
add_edge!(graph::ModelGraph,modelnodes::Vector{ModelNode}) = add_link_edge!(graph,modelnodes)

getedges(graph::ModelGraph) = graph.linkedges
getedge(graph::ModelGraph,index::Int64) = graph.linkedges[index]

function all_edges(graph::ModelGraph)
    edges = getedges(graph)
    for subgraph in graph.subgraphs
        edges = [edges;all_edges(subgraph)]
        #append!(edges,all_edges(subgraph))
    end
    return edges
end
getedge(graph::ModelGraph,nodes::Set{ModelNode}) = graph.linkedge_map[nodes]

function getedge(graph::ModelGraph,nodes::ModelNode...)
    s = Set(collect(nodes))
    return getlinkedge(graph,s)
end

function Base.getindex(graph::ModelGraph,linkedge::LinkEdge)
    return graph.edge_idx_map[linkedge]
end

# function getedge(graph::ModelGraph,vertices::Int...)
#     s = Set(collect(vertices))
#     return getlinkedge(graph,s)
# end

########################################################
# Model Management
########################################################
is_master_model(model::JuMP.Model) = haskey(model.ext,:modelgraph)
has_objective(graph::AbstractModelGraph) = graph.objective_function != zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef})
has_NLobjective(graph::AbstractModelGraph) = graph.nlp_data != nothing && graph.nlp_data.nlobj != nothing
has_subgraphs(graph::AbstractModelGraph) = !(isempty(graph.subgraphs))
has_NLlinkconstraints(graph::AbstractModelGraph) = graph.nlp_data != nothing && !(isempty(graph.nlp_data.nlconstr))
num_linkconstraints(graph::AbstractModelGraph) = length(graph.linkeqconstraints) + length(graph.linkineqconstraints)
num_linkvariables(graph::AbstractModelGraph) = length(graph.linkvariables)
num_NLlinkconstraints(graph::AbstractModelGraph) = graph.nlp_data == nothing ? 0 : length(graph.nlp_data.nlconstr)
getnumnodes(graph::AbstractModelGraph) = length(getnodes(graph))

getmasternode(graph::ModelGraph) = graph.masternode
getlinkvariables(graph::ModelGraph) = collect(values(graph.linkvariables))
getlinkconstraints(graph::ModelGraph) = collect(values(graph.linkconstraints))

#Go through subgraphs and get all linkconstraints
function all_linkconstraints(graph::AbstractModelGraph)
    links = []
    for subgraph in all_subgraphs(graph)
        append!(links,getlinkconstraints(subgraph))
    end
    append!(links,getlinkconstraints(graph))
    return links
end

function JuMP.num_variables(graph::AbstractModelGraph)
    n_master_variables = 0
    n_master_variables += JuMP.num_variables(getmasternode(graph))
    for subgraph in all_subgraphs(graph)
        n_master_variables += JuMP.num_variables(getmasternode(subgraph))
    end
    n_node_variables = sum(JuMP.num_variables.(all_nodes(graph)))
    return n_master_variables + n_node_variables
end

#JuMP Model Extenstion
####################################
# Objective
###################################
JuMP.objective_function(graph::AbstractModelGraph) = graph.objective_function
JuMP.set_objective_function(graph::AbstractModelGraph, x::JuMP.VariableRef) = JuMP.set_objective_function(graph, convert(AffExpr,x))
JuMP.set_objective_function(graph::AbstractModelGraph, func::JuMP.AbstractJuMPScalar) = graph.objective_function = func  #JuMP.set_objective_function(graph, func)

function JuMP.set_objective(graph::AbstractModelGraph, sense::MOI.OptimizationSense, func::JuMP.AbstractJuMPScalar)
    graph.objective_sense = sense
    graph.objective_function = func
end

function JuMP.objective_value(graph::AbstractModelGraph)
    objective = JuMP.objective_function(graph)
    return nodevalue(objective)
end

# nodevalue(lvref::LinkVariableRef) = JuMP.
JuMP.object_dictionary(m::ModelGraph) = m.obj_dict
JuMP.objective_sense(m::ModelGraph) = m.objective_sense

#####################################################
#  Link Variables
#  A link variable belongs to the master node on a graph.  A link variable is available for
#  model nodes to use in their constraints which 'links' them together
#####################################################
#A LinkVariableRef wraps a JuMP VariableRef.  This makes it easy to use all of the JuMP's functionality and extend it for providing link information
struct LinkVariableRef <: AbstractLinkVariableRef  #NOTE AbstractVariableRef is an AbstractJuMPScalar
    vref::JuMP.VariableRef
    idx::Int64
end

Base.broadcastable(v::LinkVariableRef) = Ref(v)
Base.copy(v::LinkVariableRef) = v
Base.:(==)(v::LinkVariableRef, w::LinkVariableRef) = v.vref.model === w.vref.model && v.vref.index == w.vref.index
JuMP.owner_model(v::LinkVariableRef) = v.vref.model
JuMP.isequal_canonical(v::LinkVariableRef, w::LinkVariableRef) = v.vref == w.vref

JuMP.variable_type(::ModelGraph) = LinkVariableRef
JuMP.set_name(v::LinkVariableRef, s::String) = JuMP.set_name(v.vref,s)
JuMP.name(v::LinkVariableRef) =  JuMP.name(v.vref)

#Add a link variable to a ModelGraph.  We need to wrap the variable in our own LinkVariableRef to work with it in constraints
function JuMP.add_variable(graph::ModelGraph, v::JuMP.AbstractVariable, name::String="")
    graph.linkvariable_index += 1

    node_vref = JuMP.add_variable(graph.masternode,v,name)  #add the variable to the master model
    link_vref = LinkVariableRef(node_vref,graph.linkvariable_index) #use underlying model variable

    graph.linkvariables[link_vref.idx] = link_vref
    JuMP.set_name(link_vref, name)

    graph.child_linkvariable_map[link_vref] = Vector{JuMP.AbstractVariableRef}()
    return link_vref
end

#Delete a link variable
function MOI.delete!(graph::ModelGraph, vref::LinkVariableRef)
    delete!(m.linkvariables, vref.idx)
    delete!(m.linkvarnames, vref.idx)
end
MOI.is_valid(graph::ModelGraph, vref::LinkVariableRef) = vref.idx in keys(graph.linkvariables)

#Link master variable to model node
#TODO: Make a linkvariable reference its master node
function link_variables!(graph::ModelGraph,lvref::LinkVariableRef,vref::JuMP.VariableRef)
    #graph = lvref.graph
    if !(vref in graph.child_linkvariable_map[lvref])
        push!(graph.child_linkvariable_map[lvref],vref)
    end
    node = getnode(vref)
    node.parent_linkvariable_map[vref] = lvref
    return nothing
end
link_variables!(graph::ModelGraph,vref::JuMP.VariableRef,lvref::LinkVariableRef) = link_variables!(graph,lvref,vref)
#link_variables!(lvref::LinkVariableRef,nvref::NodeVariableRef) = link_variables!(lvref,nvref.vref)

######################################################
# Master Constraints
######################################################
const MasterAffExpr = Union{LinkVariableRef,JuMP.GenericAffExpr{Float64,LinkVariableRef}}

struct ScalarMasterConstraint{F <: MasterAffExpr,S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::F
    set::S
    function ScalarMasterConstraint(F::MasterAffExpr,S::MOI.AbstractScalarSet)
        con = new{typeof(F),typeof(S)}(F,S)
    end
end

function JuMP.build_constraint(_error::Function,func::MasterAffExpr,set::MOI.AbstractScalarSet)
    constraint = ScalarMasterConstraint(func, set)
    return constraint
end

function JuMP.ScalarConstraint(con::ScalarMasterConstraint)
    terms = con.func.terms
    new_terms = OrderedDict([(linkvar_ref.vref,coeff) for (linkvar_ref,coeff) in terms])
    new_func = JuMP.GenericAffExpr{Float64,JuMP.VariableRef}()
    new_func.terms = new_terms
    new_func.constant = con.func.constant
    return JuMP.ScalarConstraint(new_func,con.set)
end

#Add a Master Constraint
function JuMP.add_constraint(graph::ModelGraph, con::ScalarMasterConstraint, name::String="")
    scalar_con = JuMP.ScalarConstraint(con)
    cref = JuMP.add_constraint(getmasternode(graph),scalar_con,name)          #also add to master model
    return cref
end

# Model Extras
JuMP.show_constraints_summary(::IOContext,m::ModelGraph) = ""
JuMP.show_backend_summary(::IOContext,m::ModelGraph) = ""


#####################################################
#  Link Constraints
#  A linear constraint between JuMP Models (nodes).  Link constraints can be equality or inequality.
#####################################################
struct LinkConstraintRef <: AbstractLinkConstraintRef
    graph::ModelGraph               #`model` owning the constraint
    idx::Int                        # index in `model.linkconstraints`
    linkedge::LinkEdge
end

struct LinkConstraint{F <: JuMP.AbstractJuMPScalar,S <: MOI.AbstractScalarSet} <: AbstractLinkConstraint
    func::F
    set::S
end
LinkConstraint(ref::LinkConstraintRef) = JuMP.owner_model(ref).linkconstraints[ref.idx]
LinkConstraint(con::JuMP.ScalarConstraint) = LinkConstraint(con.func,con.set)

getnodes(con::LinkConstraint) = [getnode(var) for var in keys(con.func.terms)]  #TODO: Check uniqueness.  It should be unique now that JuMP uses an OrderedDict to store terms.
getnodes(cref::LinkConstraintRef) = getnodes(cref.linkedge)
getnumnodes(con::LinkConstraint) = length(getnodes(con))


#TODO: Use macros to do this

#Add a LinkConstraint to a ModelGraph and update its LinkEdges
function add_link_equality_constraint(graph::ModelGraph,con::JuMP.ScalarConstraint,name::String = "")
    @assert isa(con.set,MOI.EqualTo)  #EQUALITY CONSTRAINTS

    graph.linkeqconstraint_index += 1
    graph.linkconstraint_index += 1

    link_con = LinkConstraint(con)    #Convert ScalarConstraint to a LinkConstraint
    modelnodes = getnodes(link_con)
    linkedge = add_link_edge!(graph,modelnodes)

    cref = LinkConstraintRef(graph, graph.linkconstraint_index,linkedge)
    JuMP.set_name(cref, name)

    push!(linkedge.linkconstraints,cref)
    graph.linkconstraints[cref.idx] = link_con
    eq_idx = graph.linkeqconstraint_index
    graph.linkeqconstraints[eq_idx] = link_con

    #Add partial linkconstraint to nodes
    for (var,coeff) in link_con.func.terms
      node = getnode(var)
      _add_to_partial_linkeqconstraint!(node,var,coeff,link_con.func.constant,link_con.set,eq_idx)
    end

    return cref
end

function add_link_inequality_constraint(graph::ModelGraph,con::JuMP.ScalarConstraint,name::String = "")
    @assert typeof(con.set) in [MOI.Interval{Float64},MOI.LessThan{Float64},MOI.GreaterThan{Float64}]

    graph.linkineqconstraint_index += 1
    graph.linkconstraint_index += 1

    link_con = LinkConstraint(con)    #Convert ScalarConstraint to a LinkConstraint
    modelnodes = getnodes(link_con)
    linkedge = add_link_edge!(graph,modelnodes)

    cref = LinkConstraintRef(graph, graph.linkconstraint_index,linkedge)
    JuMP.set_name(cref, name)

    push!(linkedge.linkconstraints,cref)
    graph.linkconstraints[cref.idx] = link_con
    ineq_idx = graph.linkineqconstraint_index
    graph.linkineqconstraints[ineq_idx] = link_con

    #Add partial linkconstraint to nodes
    for (var,coeff) in link_con.func.terms
      node = getnode(var)
      _add_to_partial_linkineqconstraint!(node,var,coeff,link_con.func.constant,link_con.set,ineq_idx)
    end

    return cref
end

function JuMP.add_constraint(graph::ModelGraph, con::JuMP.ScalarConstraint, name::String="")
    if isa(con.set,MOI.EqualTo)
        cref = add_link_equality_constraint(graph,con,name)
    else
        cref = add_link_inequality_constraint(graph,con,name)
    end
    return cref
end

#Add to a partial linkconstraint on a modelnode
function _add_to_partial_linkeqconstraint!(node::ModelNode,var::JuMP.VariableRef,coeff::Number,constant::Float64,set::MOI.AbstractScalarSet,index::Int64)
    @assert getnode(var) == node
    if haskey(node.partial_linkeqconstraints,index)
        linkcon = node.partial_linkeqconstraints[index]
        JuMP.add_to_expression!(linkcon.func,coeff,var)
    else
        new_func = JuMP.GenericAffExpr{Float64,JuMP.VariableRef}()
        new_func.terms[var] = coeff
        new_func.constant = constant
        linkcon = LinkConstraint(new_func,set)
        node.partial_linkeqconstraints[index] = linkcon
    end
end

#Add to a partial linkconstraint on a modelnode
function _add_to_partial_linkineqconstraint!(node::ModelNode,var::JuMP.VariableRef,coeff::Number,constant::Float64,set::MOI.AbstractScalarSet,index::Int64)
    @assert getnode(var) == node
    if haskey(node.partial_linkineqconstraints,index)
        linkcon = node.partial_linkineqconstraints[index]
        JuMP.add_to_expression!(linkcon.func,coeff,var)
    else
        new_func = JuMP.GenericAffExpr{Float64,JuMP.VariableRef}()
        new_func.terms[var] = coeff
        new_func.constant = constant
        linkcon = LinkConstraint(new_func,set)
        node.partial_linkineqconstraints[index] = linkcon
    end
end

function JuMP.add_constraint(graph::ModelGraph, con::JuMP.AbstractConstraint, name::String="")
    error("Cannot add constraint $con. A ModelGraph currently only supports Scalar LinkConstraints")
end

JuMP.owner_model(cref::LinkConstraintRef) = cref.graph
JuMP.constraint_type(::ModelGraph) = LinkConstraintRef
JuMP.jump_function(constraint::LinkConstraint) = constraint.func
JuMP.moi_set(constraint::LinkConstraint) = constraint.set
JuMP.shape(::LinkConstraint) = JuMP.ScalarShape()
function JuMP.constraint_object(cref::LinkConstraintRef, F::Type, S::Type)
   con = cref.graph.linkconstraints[cref.idx]
   con.func::F
   con.set::S
   return con
end
JuMP.set_name(cref::LinkConstraintRef, s::String) = JuMP.owner_model(cref).linkconstraint_names[cref.idx] = s
JuMP.name(con::LinkConstraintRef) =  JuMP.owner_model(con).linkconstraint_names[con.idx]

function MOI.delete!(graph::ModelGraph, cref::LinkConstraintRef)
    delete!(graph.linkconstraints, cref.idx)
    delete!(graph.linkconstraint_names, cref.idx)
end
MOI.is_valid(graph::ModelGraph, cref::LinkConstraintRef) = cref.idx in keys(graph.linkconstraints)

#######################################################
# HIERARCHICAL CONSTRAINTS
#######################################################
#TODO: Handle expression between link variable and node (JuMP) variables


#################################
# Optimizer
#################################
set_optimizer(graph::ModelGraph,optimizer::Union{JuMP.OptimizerFactory,AbstractGraphOptimizer,Nothing}) = graph.optimizer = optimizer


####################################
#Print Functions
####################################
function string(graph::ModelGraph)
    """
    Model Graph:
    local nodes: $(getnumnodes(graph)), total nodes: $(length(all_nodes(graph)))
    link variables: $(num_linkvariables(graph))
    local link constraints: $(num_linkconstraints(graph)), total link constraints $(length(all_linkconstraints(graph)))
    """
end
print(io::IO, graph::AbstractModelGraph) = print(io, string(graph))
show(io::IO,graph::AbstractModelGraph) = print(io,graph)
