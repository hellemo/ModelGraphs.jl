
#import JuMP:AbstractModel,setobjective,getobjectivevalue,setsolver,getvalue
##############################################################################
# ModelGraph
##############################################################################
"""
ModelGraph()

The ModelGraph Type.  Represents a graph containing models (nodes) and the linkconstraints (edges) between them.
A ModelGraph wraps a BasePlasmoGraph and can use its methods.  A ModelGraph also wraps a LinkModel object which extends a JuMP AbstractModel to provide model management functions.

"""
mutable struct ModelGraph <: AbstractModelGraph
    hypergraph::StructureGraphs.StructureGraph           #Model graph structure.  edges in the graph have references to constraints.  The graph expresses the structure of the link model
    linkmodel::LinkModel                                 #Using composition to represent a graph as a "Model".  Someday I will figure out how to do multiple inheritance.
    jump_model::Union{JuMP.AbstractModel,Nothing}        #Cache the internal serial model for the graph.  Returned if requested by the solve
end

ModelGraph() = ModelGraph(StructureGraphs.StructuredHyperGraph(),LinkModel(),nothing)
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


#################################
# Solver setters and getters
#################################
"""
setsolver(model::AbstractModelGraph,solver::AbstractMathProgSolver)

Set the graph solver to use an AbstractMathProg compliant solver
"""
#setsolver(model::AbstractModelGraph,solver::AbstractMathProgSolver) = model.linkmodel.solver = solver
set_optimizer(model::AbstractModelGraph,optimizer_factory::JuMP.OptimizerFactory) = set_optimizer(model.linkmodel,optimizer_factory)

"""
setsolver(model::AbstractModelGraph,solver::AbstractPlasmoSolver)

Set the graph solver to use an AbstractMathProg compliant solver
"""
set_optimizer(model::AbstractModelGraph,graph_solver::AbstractGraphSolver) = set_optimizer(model.linkmodel,graph_solver)

#getlinkconstraints(graph::AbstractModelGraph) = getlinkconstraints(getlinkmodel(graph))

"Get the ModelGraph solver"
#TODO Use OptimizerFactory
#getsolver(model::AbstractModelGraph) = model.linkmodel.solver

#TODO Copy a ModelGraph
# function copy(graph::AbstractModelGraph)
#     nodes = getnodes(graph)
#     edges = getedges(graph)
#     copy_graph(graph)
#     Fill in other data
# end

####################################
#Print Functions
####################################
function string(graph::AbstractModelGraph)
    "Model Graph:"*string(getlabel(graph))*"
nodes:"*string((length(getnodes(graph))))*"
link constraints (edges):"*string((getnumlinkconstraints(graph)))

end
print(io::IO, graph::AbstractModelGraph) = print(io, string(graph))
show(io::IO,graph::AbstractModelGraph) = print(io,graph)


# simple links:"*string(length(getsimplelinkconstraints(graph)))*"\n
# hyper links: "*string(length(gethyperlinkconstraints(graph)))

# #Write total objective functions for a model graph
# _setobjectivevalue(graph::AbstractModelGraph,value::Number) = graph.linkmodel.objval = value
