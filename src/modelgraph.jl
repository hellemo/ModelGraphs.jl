##############################################################################
# ModelGraph
##############################################################################
"""
ModelGraph()

The ModelGraph Type.  Represents a graph containing models (nodes) and the linkconstraints (edges) between them.
A ModelGraph wraps a BasePlasmoGraph and can use its methods.  A ModelGraph also wraps a LinkModel object which extends a JuMP AbstractModel to provide model management functions.

"""
mutable struct ModelGraph <: AbstractModelGraph
    hypergraph::StructureGraphs.StructureGraph           #Model graph structure.  Edges and Nodes in the graph have references to LinkConstraints.  The graph expresses the structure of the link model
    linkmodel::AbstractLinkModel                         #A ModelGraph and its underlying LinkModel maintain references to eachother.
    jump_model::Union{JuMP.AbstractModel,Nothing}        #Cache the internal serial model for the graph if using an MOI solver.  Returned if requested by the solve.
    function ModelGraph()
        graph = new()
        graph.hypergraph = StructureGraphs.StructuredHyperGraph()
        graph.linkmodel = LinkModel(graph)
        graph.jump_model = nothing
        return graph
    end
end

#ModelGraph() = ModelGraph(StructureGraphs.StructuredHyperGraph(),LinkModel(),nothing)
StructureGraphs.getstructuregraph(graph::ModelGraph) = graph.hypergraph

getlinkmodel(graph::AbstractModelGraph) = graph.linkmodel
getnumlinkconstraints(graph::AbstractModelGraph) = length(graph.linkmodel.linkconstraints)
hasobjective(graph::AbstractModelGraph) = getlinkmodel(graph).objective_function != zero(JuMP.GenericAffExpr{Float64, JuMP.AbstractVariableRef})

"Set the objective of a ModelGraph"
#TODO. Write objective methods for the LinkModel
JuMP.set_objective_function(graph::AbstractModelGraph, sense::MOI.OptimizationSense, x::JuMP.VariableRef) = JuMP.set_objective_function(graph.linkmodel, sense, convert(AffExpr,x))

"Get the ModelGraph objective value"
JuMP.getobjectivevalue(graph::AbstractModelGraph) = getobjectivevalue(graph.linkmodel)

"Get the current created JuMP model for the ModelGraph.  Only created when solving using a JuMP compliant solver."
getinternaljumpmodel(graph::AbstractModelGraph) = graph.serial_model

getgraphvariables(graph::AbstractModelGraph) = getgraphvariables(getlinkmodel(graph))
getgraphconstraints(graph::AbstractModelGraph) = getgraphconstraints(getlinkmodel(graph))

"""
get_all_linkconstraints(graph::AbstractModelGraph)

Get a list containing every link constraint in the graph, including its subgraphs
"""
function get_all_linkconstraints(graph::AbstractModelGraph)
    links = []
    for subgraph in getsubgraphlist(graph)
        append!(links,getlinkconstraints(subgraph))
    end
    append!(links,getlinkconstraints(graph))
    return links
end

#################################
# Solver setters and getters
#################################
#TODO Figure out how this works in JuMP 0.19.  Might have to hook into MOI here.
#set_optimizer(model::AbstractModelGraph,graph_solver::AbstractGraphSolver) = set_optimizer(model.linkmodel,graph_solver)

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
