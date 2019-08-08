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
    hypergraph::NestedHyperGraph

    #Store master variables and constraints in a stand-alone JuMP Model
    mastermodel::JuMP.Model

    #Map from hypernodes and hyperedges to model nodes and link edges
    modelnodes::Dict{HyperNode,ModelNode}
    linkedges::Dict{HyperEdge,LinkEdge}

    #Link variables
    linkvariables::Dict{Int,JuMP.AbstractVariable}                        #Link Variable reference.  These variables are also in the mastermodel.
    linkvariablemap::Dict{JuMP.AbstractVariable,JuMP.AbstractVariable}    #map of link variables in master model to corresponding variables in ModelNodes.
    linkvariablenames::Dict{Int,String}

    #Link constraints
    linkconstraints::Dict{Int,AbstractLinkConstraint}                     #Link constraint.  Defined over variables in ModelNodes.
    linkconstraintnames::Dict{Int,String}

    #Objective
    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar

    #Optimizer
    optimizer::Union{JuMP.OptimizerFactory,AbstractGraphOptimizer,Nothing}

    link_variable_index::Int
    link_constraint_index::Int                              #keep track of master and link constraints

    obj_dict::Dict{Symbol,Any}

    #Nonlinear Link Constraints
    nlp_data::JuMP._NLPData

    #Constructor
    function ModelGraph()
        model = new(NestedHyperGraph(),
                    JuMP.Model(),
                    Dict{HyperNode,ModelNode}(),
                    Dict{HyperEdge,LinkEdge}(),
                    Dict{Int, JuMP.AbstractVariable}(),
                    Dict{JuMP.AbstractVariable, JuMP.AbstractVariable}(),
                    Dict{Int, JuMP.AbstractConstraint}(),
                    Dict{Int, AbstractLinkConstraint}(),
                    Dict{Int, String}(),
                    Dict{Int, String}(),
                    Dict{Int, String}(),
                    MOI.FEASIBILITY_SENSE,
                    zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef}),
                    #nothing,
                    0,
                    0,
                    Dict{Symbol,Any}()
                    )
        return model
    end
end

########################################################
# HyperGraph Interface
########################################################
gethypergraph(modelgraph::ModelGraph) = modelgraph.hypergraph
function add_subgraph!(graph::ModelGraph,subgraph::ModelGraph)
    hypergraph = gethypergraph(graph)
    sub_hypergraph = gethypergraph(subgraph)
    add_subgraph!(hypergraph,subgraph)
    return nothing
end

getmodelnode(graph::ModelGraph,node::HyperNode) = graph.modelnodes[hypernode]

function getnode(graph::ModelGraph,index::Int64)
    hypernode = getnode(gethypergraph(graph),index)
    return getmodelnode(graph,hypernode)
end

function getnodes(graph::ModelGraph)
    return getmodelnode.(getnodes(graph.hypergraph))
end

function getindex(graph::ModelGraph,node::ModelNode)
    hypernode = gethypernode(node)
    return getindex(gethypergraph(graph),hypernode)
end


getlinkedge(graph::ModelGraph,hyperedge::HyperEdge) = graph.linkedges[hyperedge]

function getlinkedge(graph::ModelGraph,index::Int64)
    hyperedge = gethyperedge(gethypergraph(graph),index)
    return getlinkedge(graph,hyperedge)
end

function getlinkedge(graph::ModelGraph,vertices::Int...)
    hyperedge = gethyperedge(graph.hypergraph,vertices...)
    return getlinkedge(graph,hyperedge)
end

function getlinkedges(graph::ModelGraph)
    return getlinkedge.(gethyperedges(graph.hypergraph))
end

function getindex(graph::ModelGraph,linkedge::LinkEdge)
    hyperedge = gethyperedge(linkedge)
    return getindex(gethypergraph(graph),hyperedge)
end
########################################################
# Model Management
########################################################

has_objective(graph::AbstractModelGraph) = graph.objective_function != zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef})
get_num_linkconstraints(graph::AbstractModelGraph) = length(graph.linkconstraints)
get_num_NLlinkconstraints(graph::AbstractModelGraph) = graph.nlp_data == nothing? 0 : length(graph.nlp_data.nlconstr)

getmastermodel(graph) = graph.mastermodel

get_linkvariables(m::ModelGraph) = collect(values(graph.link_variables))
get_linkconstraints(graph::ModelGraph) = collect(values(model.linkconstraints))

function get_all_linkconstraints(graph::AbstractModelGraph)
    links = []
    for subgraph in subgraphs(graph)
        append!(links,getlinkconstraints(subgraph))
    end
    append!(links,getlinkconstraints(graph))
    return links
end

#JuMP Model Extenstion
#TODO. Write more objective methods
JuMP.set_objective_function(graph::AbstractModelGraph, sense::MOI.OptimizationSense, x::JuMP.VariableRef) = JuMP.set_objective_function(graph, sense, convert(AffExpr,x))
JuMP.objective_value(graph::AbstractModelGraph) = getobjectivevalue(graph.linkmodel)

JuMP.object_dictionary(m::ModelGraph) = m.objdict
JuMP.objective_sense(m::ModelGraph) = m.objective_sense

#####################################################
#  Link Variables
#  A link variable belongs to the master node on a graph.  A link variable is available for
#  model nodes to use in their constraints which 'links' them together
#####################################################
#Reference to a Link Variable.  These variables can be referenced in constraints from JuMP models
struct LinkVariableRef <: JuMP.AbstractVariableRef  #NOTE AbstractVariableRef is an AbstractJuMPScalar
    model::JuMP.Model
    idx::Int                    #index in the master model
end

Base.broadcastable(v::LinkVariableRef) = Ref(v)
Base.copy(v::LinkVariableRef) = v
Base.:(==)(v::LinkVariableRef, w::LinkVariableRef) = v.model === w.model && v.idx == w.idx
JuMP.owner_model(v::LinkVariableRef) = v.model
JuMP.isequal_canonical(v::LinkVariableRef, w::LinkVariableRef) = v == w

JuMP.variable_type(::ModelGraph) = LinkVariableRef

JuMP.set_name(v::LinkVariableRef, s::String) = JuMP.owner_model(v).linkvarnames[v.idx] = s
JuMP.name(v::LinkVariableRef) =  JuMP.owner_model(v).linkvarnames[v.idx]

#Add a link variable to a ModelGraph.  We need to wrap the variable in our own LinkVariableRef to work with it in constraints
function JuMP.add_variable(graph::ModelGraph, v::JuMP.AbstractVariable, name::String="")
    graph.link_var_index += 1
    m = getmastermodel(graph)
    vref = LinkVariableRef(m, graph.link_var_index)
    JuMP.add_variable(graph.mastermodel,v,name)  #add the variable to the master model
    graph.linkvariables[vref.idx] = v
    JuMP.set_name(vref, name)
    graph.linkvariablemap[vref] = Vector{JuMP.AbstractVariable}()
    return vref
end

#Delete a link variable
function MOI.delete!(graph::ModelGraph, vref::LinkVariableRef)
    delete!(m.linkvariables, vref.idx)
    delete!(m.linkvarnames, vref.idx)
end
MOI.is_valid(graph::ModelGraph, vref::LinkVariableRef) = vref.idx in keys(graph.linkvariables)
JuMP.num_link_variables(graph::ModelGraph) = length(graph.linkvariables)

#Link master variable to model node
function link_variables!(lvref::LinkVariableRef,vref::JuMP.VariableRef)
    graph = lvref.modelgraph
    if !(vref in graph.linkvariablemap[lvref])
        push!(graph.linkvariablemap[lvref],vref)
    end
    node = getnode(vref)
    node.linkvariablemap[vref] = lvref
    return nothing
end

######################################################
# Master Constraints
######################################################
const MasterAffExpr = JuMP.GenericAffExpr{Float64,LinkVariableRef}
struct ScalarMasterConstraint{F <: MasterAffExpr,S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::F
    set::S
    function ScalarMasterConstraint(F::AbstractJuMPScalar,S::MOI.AbstractScalarSet)
        con = new{typeof(F),typeof(S)}(F,S)
    end
end

function JuMP.build_constraint(_error::Function,func::MasterAffExpr,set::MOI.AbstractScalarSet)
    constraint = ScalarMasterConstraint(func, set)
    return constraint
end

#Add a Master Constraint
function JuMP.add_constraint(graph::ModelGraph, con::ScalarLinkVarConstraint, name::String="")
    #graph.linkvar_constraint_index += 1
    #cref = LinkVarConstraintRef(m, m.linkvar_constraint_index)
    cref = JuMP.add_constraint(getmastermodel(graph),con,name)          #also add to master model
    return cref
end

# Model Extras
JuMP.show_constraints_summary(::IOContext,m::ModelGraph) = ""
JuMP.show_backend_summary(::IOContext,m::ModelGraph) = ""


#####################################################
#  Link Constraints
#  A linear constraint between JuMP Models (nodes).  Link constraints can be equality or inequality.
#####################################################
struct LinkConstraintRef
    graph::ModelGraph               #`model` owning the constraint
    idx::Int                        # index in `model.linkconstraints`
    linkedge::LinkEdge
end

JuMP.constraint_type(::ModelGraph) = LinkConstraintRef #GraphConstraintRef
JuMP.owner_model(con::LinkConstraintRef) = con.graph

struct LinkConstraint{F <: JuMP.AbstractJuMPScalar,S <: MOI.AbstractScalarSet} <: AbstractLinkConstraint
    func::F
    set::S
end
LinkConstraint(ref::LinkConstraintRef) = JuMP.owner_model(ref).linkconstraints[ref.idx]
LinkConstraint(con::JuMP.ScalarConstraint) = LinkConstraint(con.func,con.set)

function getnodes(con::LinkConstraint)
    return [getnode(var) for var in keys(con.func.terms)]   #TODO: Check uniqueness.  It should be unique now that JuMP uses an OrderedDict to store terms.
    #return getnodes(con.link_edge)
    #return getnodes(con.hyperedge)
    #return map(n -> getnode(con.graph,n),con.node_indices)
end
getnodes(cref::LinkConstraintRef) = getnodes(cref.linkedge)
getnumnodes(con::LinkConstraint) = length(getnodes(con))

#Add a LinkConstraint to a ModelGraph and Update its LinkEdges
function JuMP.add_constraint(graph::ModelGraph, con::JuMP.ScalarConstraint, name::String="")
    graph.link_constraint_index += 1
    link_con = LinkConstraint(con)      #convert ScalarConstraint to a LinkConstraint

    #Setup graph information
    hypergraph = gethypergraph(graph)
    hypernode_indices = sort(unique([getindex(hypergraph,getnode(var)) for var in keys(con.func.terms)]))
    model_nodes = getnode.(graph,hypernode_indices)

    #hyperedge = getedge(hypergraph,hypernode_indices)
    link_edge = add_link_edge!(graph,model_nodes)

    cref = LinkConstraintRef(graph, graph.link_constraint_index,link_edge)
    graph.linkconstraints[cref.idx] = link_con

    #graph.linkconstraint_hyperedge_map[link_con] = hyperedge
    JuMP.set_name(cref, name)

    #Add LinkEdges to the Graph
    #add_link_edge!(graph,cref)

    return cref
end

function JuMP.add_constraint(graph::ModelGraph, con::JuMP.AbstractConstraint, name::String="")
    error("A ModelGraph only supports LinkConstraints")
end

jump_function(constraint::LinkConstraint) = constraint.func
moi_set(constraint::LinkConstraint) = constraint.set
shape(::LinkConstraint) = JuMP.ScalarShape()


#################################
# Optimizer
#################################
set_optimizer(graph::AbstractModelGraph,optimizer::Union{JuMP.OptimizerFactory,AbstractGraphOptimizer,Nothing}) = graph.optimizer = optimizer


####################################
#Print Functions
####################################
function string(graph::ModelGraph)
    "Model Graph:"*string(getlabel(graph))*"
nodes:"*string((length(getnodes(graph))))*"
link constraints (edges):"*string((getnumlinkconstraints(graph)))
end
print(io::IO, graph::AbstractModelGraph) = print(io, string(graph))
show(io::IO,graph::AbstractModelGraph) = print(io,graph)




# #Caching data on the LinkConstraint
# graph::AbstractModelGraph
# node_indices::Vector{Int64}   #Should be easier to just reference a hyperedge
# hyperedge::HyperEdge  #reference to a hyperedge



# function LinkConstraint(con::JuMP.ScalarConstraint,graph::AbstractModelGraph)
#     node_indices = sort(unique([getindex(graph,getnode(var)) for var in keys(con.func.terms)]))
#     return LinkConstraint(con.func,con.set,graph,node_indices)
# end

#Constraint Reference to a Constraint containing LinkVariables (i.e. a mastermodel constraint)
# struct LinkVarConstraintRef <: AbstractGraphConstraintRef
#     graph::ModelGraph               #`model` owning the constraint
#     idx::Int                        #index in `model.graphconstraints`
# end
# JuMP.constraint_type(::ModelGraph) = AbstractGraphConstraintRef #GraphConstraintRef
# JuMP.owner_model(con::LinkVarConstraintRef) = con.graph

# function JuMP.constraint_object(cref::LinkVarConstraintRef, F::Type, S::Type)
#    con = cref.graph.allconstraints[cref.idx]
#    con.func::F
#    con.set::S
#    return con
# end

# function MOI.delete!(graph::ModelGraph, cref::LinkVarConstraintRef)
#     # if typeof(constraint_object(cref)) == LinkConstraint
#     #     delete!(m.linkconstraints, cref.idx)
#     #     delete!(m.linkconnames, cref.idx)
#     #elseif typeof(constraint_object(cref)) == JuMP.AbstractConstraint
#     delete!(graph.mastermodel, cref)
#     #else
#         #error("Constraint reference refers to constraint type that was not recognized")
#     #end
# end
#
# MOI.is_valid(m::ModelGraph, cref::LinkVarConstraintRef) = cref.idx in keys(m.)#keys(m.graphconstraints) || cref.idx in keys(m.linkconstraints)
#
# JuMP.set_name(con::LinkVarConstraintRef, s::String) = JuMP.owner_model(con).graphconstraintnames[con.idx] = s
# JuMP.name(con::GraphConstraintRef) =  JuMP.owner_model(con).graphconstraintnames[con.idx]
