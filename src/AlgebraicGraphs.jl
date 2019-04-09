module AlgebraicGraphs

using StructureGraphs
#import LightGraphs

using Requires
using Distributed
using LinearAlgebra

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

import JuMP
import JuMP: AbstractModel, AbstractConstraint, AbstractJuMPScalar, ConstraintRef
import Base: ==,show,print,string,getindex,copy

#Model Graph Constructs
export AbstractModelGraph, ModelGraph, SolutionGraph, JuMPGraph,

#Graph transformations
ModelBipartiteGraph, NodeUnipartiteGraph, LinkUnipartiteGraph,

#Nodes and Edges,
ModelNode,LinkingEdge,

#LinkModel Types
LinkConstraint, GraphConstraint, GraphVariableRef, GraphConstraintRef,

#Solvers
AbstractGraphSolver,BendersSolver,LagrangeSolver,

#re-export base functions
addnode!,add_node!,getnodes,getedges,collectnodes,

#Model functions
setmodel,setsolver,setmodel!,resetmodel,is_nodevar,getmodel,getsolver,hasmodel,
getnumnodes, getobjectivevalue, getinternalgraphmodel,

#Link Constraints
addlinkconstraint, getlinkreferences, getlinkconstraints, getsimplelinkconstraints, gethyperlinkconstraints, get_all_linkconstraints,


#Graph Transformation functions
aggregate!,create_aggregate_model,create_partitioned_model_graph,create_lifted_model_graph,getbipartitegraph,getunipartitegraph,partition,label_propagation,

#JuMP Interface functions
create_jump_graph_model,
getgraph,getnodevariables,getnodevariable,getnodevariablemap,getnodeobjective,getnodeconstraints,getnodedata,is_graphmodel,

#solve handles
solve_jump,bendersolve,solve,

#Solution management
getsolution,setsolution,setvalue,getvalue,

#macros
@linkconstraint,@graphobjective

#Abstract Types
abstract type AbstractModelGraph <: AbstractStructureGraph end
abstract type AbstractModelNode <: AbstractStructureNode end
abstract type AbstractLinkingEdge  <: AbstractStructureEdge end
abstract type AbstractGraphSolver end

include("linkmodel.jl")  #A JuMP extension model to manage special variables and constraints

include("modelgraph.jl") #The ModelGraph

include("modelnode.jl")  #ModelGraph nodes

include("linkconstraint.jl")  #LinkConstraints which create graph topology

include("modeledge.jl")   #ModelGraph edges

include("macros.jl")      #@linkconstraint, @graphobjective

include("jumpgraph.jl")   #An aggregated JuMP model

include("solve.jl")       #Aggregate and solve with an MOI Solver
#
# include("solution.jl")  #SolutionGraph
#

#
# include("aggregation.jl")  #Aggregate pieces of a ModelGraph
#
# include("community_detection.jl")
#
# include("graph_transformations/modeltree.jl")
#
# include("graph_transformations/pipstree.jl")
#
# include("graph_transformations/partite_graphs.jl")   #Transformations that facilitate graph analysis (partitioning and community detection)
#
# include("graph_transformations/graph_transformation.jl")
#
# function __init__()
#     @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("extras/plots.jl")
#     @require Metis = "2679e427-3c69-5b7f-982b-ece356f1e94b" include("extras/metis.jl")
#     @require CommunityDetection = "d427f087-d71a-5a1b-ace0-b93392eea9ff" include("extras/communities.jl")
# end


# if haskey(Pkg.installed(),"MPI")
# #External Solver Interfaces
#     include("solver_interfaces/wrapped_solvers.jl")
#
#     include("solver_interfaces/plasmoPipsNlpInterface.jl")
#     using .PlasmoPipsNlpInterface
# end

end
