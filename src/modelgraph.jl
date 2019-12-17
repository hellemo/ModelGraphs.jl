##############################################################################
# ModelGraph
##############################################################################
"""
ModelGraph()

The ModelGraph Type.  Represents a graph containing models (nodes) and the linkconstraints (edges) between them.
A ModelGraph wraps a BasePlasmoGraph and can use its methods.  A ModelGraph also wraps a LinkModel object which extends a JuMP AbstractModel to provide model management functions.

"""
mutable struct ModelGraph <: AbstractModelGraph
    #Nested HyperGraph represents the problem structure
    hypergraph::HyperGraph

    #Store master variables and constraints in a stand-alone JuMP Model
    mastermodel::JuMP.Model #TODO: Make this a master node.  It will be possible to have link constraints that connect a master node to other nodes

    #Map from hypernodes and hyperedges to model nodes and link edges
    modelnodes::Dict{HyperNode,ModelNode}
    linkedges::Dict{HyperEdge,LinkEdge}         #local link edges
    subgraphs::Vector{AbstractModelGraph}

    #Link variables
    masterlinkvariables::Dict{AbstractLinkVariableRef,AbstractLinkVariableRef}      #Link Variables from higher level master model
    linkvariables::OrderedDict{Int,AbstractLinkVariableRef}                                #Link Variables in master model
    linkvariable_map::Dict{AbstractLinkVariableRef,Vector{JuMP.VariableRef}}        #Map of link variables in master model to corresponding variables in ModelNodes.
    linkvariable_names::Dict{Int,String}

    #Link constraints
    linkconstraints::OrderedDict{Int,AbstractLinkConstraint}                     #Link constraint.  Defined over variables in ModelNodes.
    linkconstraint_names::Dict{Int,String}

    #Objective
    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar

    #Optimizer
    optimizer::Union{JuMP.OptimizerFactory,AbstractGraphOptimizer,Nothing}

    linkvariable_index::Int
    linkconstraint_index::Int                              #keep track of master and link constraints

    obj_dict::Dict{Symbol,Any}

    #Nonlinear Link Constraints
    nlp_data::Union{Nothing,JuMP._NLPData}

    #Constructor
    function ModelGraph()
        modelgraph = new(HyperGraph(),
                    JuMP.Model(),
                    Dict{HyperNode,ModelNode}(),
                    Dict{HyperEdge,LinkEdge}(),
                    Vector{AbstractModelGraph}(),
                    Dict{AbstractLinkVariableRef,AbstractLinkVariableRef}(),
                    OrderedDict{Int, JuMP.AbstractVariable}(),
                    Dict{JuMP.AbstractVariable, JuMP.AbstractVariable}(),
                    Dict{Int,String}(),
                    OrderedDict{Int, AbstractLinkConstraint}(),
                    Dict{Int, String}(),
                    MOI.FEASIBILITY_SENSE,
                    zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef}),
                    nothing,
                    0,
                    0,
                    Dict{Symbol,Any}(),
                    nothing
                    )

        modelgraph.mastermodel.ext[:modelgraph] = modelgraph

        return modelgraph
    end
end


########################################################
# HyperGraph Interface
########################################################
gethypergraph(modelgraph::ModelGraph) = modelgraph.hypergraph
function add_subgraph!(graph::ModelGraph,subgraph::ModelGraph)
    hypergraph = gethypergraph(graph)
    sub_hypergraph = gethypergraph(subgraph)
    add_subgraph!(hypergraph,sub_hypergraph)
    push!(graph.subgraphs,subgraph)

    for node in getnodes(sub_hypergraph)
        graph.modelnodes[node] = getnode(subgraph,node)
    end

    # for hyperedge in getedges(sub_hypergraph)
    #     graph.sublinkedges[hyperedge] = getlinkedge(subgraph,hyperedge)
    # end

    return nothing
end

subgraphs(modelgraph) = modelgraph.subgraphs

getmodelnode(graph::ModelGraph,hypernode::HyperNode) = graph.modelnodes[hypernode]
getnode(graph::ModelGraph,hypernode::HyperNode) = getmodelnode(graph,hypernode)

function getnode(graph::ModelGraph,index::Int64)
    hypernode = getnode(gethypergraph(graph),index)
    return getmodelnode(graph,hypernode)
end
function getnodes(graph::ModelGraph)
    return map(x -> getnode(graph,x),getnodes(graph.hypergraph))    #return ge   tmodelnode.(getnodes(graph.hypergraph))
end



function Base.getindex(graph::ModelGraph,node::ModelNode)
    hypernode = gethypernode(node)
    return getindex(gethypergraph(graph),hypernode)
end


getlocallinkedge(graph::ModelGraph,hyperedge::HyperEdge) = graph.linkedges[hyperedge]
getsublinkedge(graph::ModelGraph,hyperedge::HyperEdge) = graph.sublinkedges[hyperedge]
function getlinkedge(graph::ModelGraph,hyperedge::HyperEdge)
    return graph.linkedges[hyperedge]
end

function findlinkedge(graph::ModelGraph,hyperedge::HyperEdge)
    if haskey(graph.linkedges,hyperedge)
        return graph.linkedges[hyperedge]
    else
        for subgraph in subgraphs(graph)
            if haskey(subgraph.linkedges,hyperedge)
                return subgraph.linkedges[hyperedge]
            end
        end
    end
end

function getlinkedge(graph::ModelGraph,index::Int64)
    hyperedge = gethyperedge(gethypergraph(graph),index)
    return getlinkedge(graph,hyperedge)
end

function getlinkedge(graph::ModelGraph,vertices::Int...)
    hyperedge = gethyperedge(graph.hypergraph,vertices...)
    return getlinkedge(graph,hyperedge)
end

function getlinkedges(graph::ModelGraph)
    return map(x -> getlinkedge(graph,x),gethyperedges(graph.hypergraph))
    #return getlinkedge.(gethyperedges(graph.hypergraph))
end

function getindex(graph::ModelGraph,linkedge::LinkEdge)
    hyperedge = gethyperedge(linkedge)
    return getindex(gethypergraph(graph),hyperedge)
end
########################################################
# Model Management
########################################################
is_master_model(model::JuMP.Model) = haskey(model.ext,:modelgraph)
has_objective(graph::AbstractModelGraph) = graph.objective_function != zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef})
has_NLobjective(graph::AbstractModelGraph) = graph.nlp_data != nothing && graph.nlp_data.nlobj != nothing
has_subgraphs(graph::AbstractModelGraph) = !(isempty(graph.subgraphs))
has_NLlinkconstraints(graph::AbstractModelGraph) = graph.nlp_data != nothing && !(isempty(graph.nlp_data.nlconstr))
num_linkconstraints(graph::AbstractModelGraph) = length(graph.linkconstraints)
num_linkvariables(graph::AbstractModelGraph) = length(graph.linkvariables)
num_NLlinkconstraints(graph::AbstractModelGraph) = graph.nlp_data == nothing ? 0 : length(graph.nlp_data.nlconstr)
getnumnodes(graph::AbstractModelGraph) = length(getnodes(gethypergraph(graph)))

getmastermodel(graph) = graph.mastermodel

#NOTE: Order is random
getlinkvariables(graph::ModelGraph) = collect(values(graph.linkvariables))
getlinkconstraints(graph::ModelGraph) = collect(values(graph.linkconstraints))

#Go through subgraphs and get all linkconstraints
function all_linkconstraints(graph::AbstractModelGraph)
    links = []
    for subgraph in subgraphs(graph)
        append!(links,getlinkconstraints(subgraph))
    end
    append!(links,getlinkconstraints(graph))
    return links
end

function JuMP.num_variables(graph::AbstractModelGraph)
    n_master_variables = 0
    n_master_variables += JuMP.num_variables(getmastermodel(graph))
    for subgraph in subgraphs(graph)
        n_master_variables += JuMP.num_variables(getmastermodel(subgraph))
    end

    n_node_variables = sum(JuMP.num_variables.(getnodes(graph)))

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
    graph::ModelGraph
    idx::Int64
end

Base.broadcastable(v::LinkVariableRef) = Ref(v)
Base.copy(v::LinkVariableRef) = v
Base.:(==)(v::LinkVariableRef, w::LinkVariableRef) = v.vref.model === w.vref.model && v.vref.index == w.vref.index
JuMP.owner_model(v::LinkVariableRef) = v.graph.mastermodel
JuMP.isequal_canonical(v::LinkVariableRef, w::LinkVariableRef) = v.vref == w.vref

JuMP.variable_type(::ModelGraph) = LinkVariableRef

JuMP.set_name(v::LinkVariableRef, s::String) = JuMP.set_name(v.vref,s)
JuMP.name(v::LinkVariableRef) =  JuMP.name(v.vref)

# JuMP.set_name(v::LinkVariableRef, s::String) = JuMP.owner_model(v).linkvarnames[v.idx] = s
# JuMP.name(v::LinkVariableRef) =  JuMP.owner_model(v).linkvarnames[v.idx]

#Add a link variable to a ModelGraph.  We need to wrap the variable in our own LinkVariableRef to work with it in constraints
function JuMP.add_variable(graph::ModelGraph, v::JuMP.AbstractVariable, name::String="")
    graph.linkvariable_index += 1

    j_vref = JuMP.add_variable(graph.mastermodel,v,name)  #add the variable to the master model
    vref = LinkVariableRef(j_vref,graph,graph.linkvariable_index)

    graph.linkvariables[j_vref.index.value] = vref
    JuMP.set_name(vref, name)

    graph.linkvariable_map[vref] = Vector{JuMP.VariableRef}()
    return vref
end

#Delete a link variable
function MOI.delete!(graph::ModelGraph, vref::LinkVariableRef)
    delete!(m.linkvariables, vref.idx)
    delete!(m.linkvarnames, vref.idx)
end
MOI.is_valid(graph::ModelGraph, vref::LinkVariableRef) = vref.idx in keys(graph.linkvariables)

#Link master variable to model node
function link_variables!(lvref::LinkVariableRef,vref::JuMP.VariableRef)
    graph = lvref.graph
    if !(vref in graph.linkvariable_map[lvref])
        push!(graph.linkvariable_map[lvref],vref)
    end
    node = getnode(vref)
    node.linkvariablemap[vref] = lvref
    return nothing
end
link_variables!(vref::JuMP.VariableRef,lvref::LinkVariableRef) = link_variables!(lvref,vref)

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
    cref = JuMP.add_constraint(getmastermodel(graph),scalar_con,name)          #also add to master model
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

function getnodes(con::LinkConstraint)
    return [getnode(var) for var in keys(con.func.terms)]   #TODO: Check uniqueness.  It should be unique now that JuMP uses an OrderedDict to store terms.
end
#return getnodes(con.link_edge)
#return getnodes(con.hyperedge)
#return map(n -> getnode(con.graph,n),con.node_indices)
getnodes(cref::LinkConstraintRef) = getnodes(cref.linkedge)
getnumnodes(con::LinkConstraint) = length(getnodes(con))

#Add a LinkConstraint to a ModelGraph and update its LinkEdges
function JuMP.add_constraint(graph::ModelGraph, con::JuMP.ScalarConstraint, name::String="")
    graph.linkconstraint_index += 1
    link_con = LinkConstraint(con)      #convert ScalarConstraint to a LinkConstraint

    #Setup graph information
    hypergraph = gethypergraph(graph)
    hypernodes = sort(unique([getindex(hypergraph,getnode(var).hypernode) for var in keys(con.func.terms)]))
    modelnodes = [getnode(graph,index) for index in hypernodes]

    linkedge = add_link_edge!(graph,modelnodes)

    cref = LinkConstraintRef(graph, graph.linkconstraint_index,linkedge)
    push!(linkedge.linkconstraints,cref)
    graph.linkconstraints[cref.idx] = link_con
    JuMP.set_name(cref, name)

    for (var,coeff) in link_con.func.terms
      node = getnode(var)
      _add_to_partial_linkconstraint!(node,var,coeff,link_con.func.constant,link_con.set,cref.idx)
    end
    return cref
end

#Add to a partial linkconstraint on a modelnode
function _add_to_partial_linkconstraint!(node::ModelNode,var::JuMP.VariableRef,coeff::Number,constant::Float64,set::MOI.AbstractScalarSet,index::Int64)
    @assert getnode(var) == node
    if haskey(node.partial_linkconstraints,index)
        linkcon = node.partial_linkconstraints[index]
        JuMP.add_to_expression!(linkcon.func, coeff,var)
    else
        new_func = JuMP.GenericAffExpr{Float64,JuMP.VariableRef}()
        new_func.terms[var] = coeff
        new_func.constant = constant
        linkcon = LinkConstraint(new_func,set)
        node.partial_linkconstraints[index] = linkcon
    end
end

function JuMP.add_constraint(graph::ModelGraph, con::JuMP.AbstractConstraint, name::String="")
    error("A ModelGraph only supports LinkConstraints")
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


#################################
# Optimizer
#################################
set_optimizer(graph::AbstractModelGraph,optimizer::Union{JuMP.OptimizerFactory,AbstractGraphOptimizer,Nothing}) = graph.optimizer = optimizer


####################################
#Print Functions
####################################
function string(graph::ModelGraph)
    """
    Model Graph:
    model nodes: $(getnumnodes(graph))
    link variables: $(num_linkvariables(graph))
    link constraints: $(num_linkconstraints(graph))
    """
end
print(io::IO, graph::AbstractModelGraph) = print(io, string(graph))
show(io::IO,graph::AbstractModelGraph) = print(io,graph)
