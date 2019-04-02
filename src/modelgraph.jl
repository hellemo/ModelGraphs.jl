import StructureGraphs:add_node!,add_edge!,create_node,create_edge,getnode,getnodes
import Base:show,print,string,getindex,copy
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
    hypergraph::StructuredHyperGraph                     #Model graph structure.  edges in the graph have references to constraints.  The graph expresses the structure of the link model
    linkmodel::LinkModel                                 #Using composition to represent a graph as a "Model".  Someday I will figure out how to do multiple inheritance.
    jump_model::Union{JuMP.AbstractModel,Nothing}        #Cache the internal serial model for the graph.  Returned if requested by the solve
end

ModelGraph() = ModelGraph(StructuredHyperGraph(),LinkModel(),nothing)

# #Write total objective functions for a model graph
# _setobjectivevalue(graph::AbstractModelGraph,value::Number) = graph.linkmodel.objval = value

"Set the objective of a ModelGraph"
JuMP.setobjective(graph::AbstractModelGraph, sense::MOI.OptimizationSense, x::JuMP.Variable) = setobjective(graph.linkmodel, sense, convert(AffExpr,x))

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

"""
setsolver(model::AbstractModelGraph,solver::AbstractPlasmoSolver)

Set the graph solver to use an AbstractMathProg compliant solver
"""
#setsolver(model::AbstractModelGraph,solver::AbstractPlasmoSolver) = model.linkmodel.solver = solver

"Get the ModelGraph solver"
#getsolver(model::AbstractModelGraph) = model.linkmodel.solver

########################################
# Link Constraints
########################################
#Store link constraint in the given graph.  Store a reference to the linking constraint on the nodes which it links

#TODO
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
    "Model Graph\ngraph_id: "*string(getlabel(graph))*"\nnodes:"*string((length(getnodes(graph))))*"\nsimple links:"*string(length(getsimplelinkconstraints(graph)))*"\nhyper links: "*string(length(gethyperlinkconstraints(graph)))
end
print(io::IO, graph::AbstractModelGraph) = print(io, string(graph))
show(io::IO,graph::AbstractModelGraph) = print(io,graph)
