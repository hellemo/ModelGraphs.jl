module AlgebraicGraphs

using StructureGraphs
import LightGraphs

using Requires
using Distributed
using LinearAlgebra
using DataStructures
using SparseArrays

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

import JuMP
import JuMP: AbstractModel, AbstractConstraint, AbstractJuMPScalar, ConstraintRef
import Base: ==,show,print,string,getindex,copy

#Model Graph Constructs
export AbstractModelGraph, ModelGraph, SolutionGraph, JuMPGraph,

#Graph Projections
ModelBipartiteGraph, NodeUnipartiteGraph, LinkUnipartiteGraph,

#Nodes and Edges,
ModelNode,LinkingEdge,

#LinkModel Types
LinkConstraint, GraphConstraint, GraphVariableRef, GraphConstraintRef,

#Solvers
AbstractGraphSolver,BendersSolver,LagrangeSolver,

#re-export base functions
addnode!,add_node!,getnode,getnodes,getedges,collectnodes,

#Model functions
setmodel,setsolver,setmodel!,resetmodel,is_nodevar,getmodel,hasmodel,
getnumnodes, getobjectivevalue, getinternalgraphmodel,

#Link Constraints
addlinkconstraint, getlinkreferences, getlinkconstraints, getsimplelinkconstraints, gethyperlinkconstraints, get_all_linkconstraints,

#Partitioning
partition,

#Aggregation
create_aggregate_model,create_aggregate_graph,

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

#ModelGraph
abstract type AbstractModelGraph <: AbstractStructureGraph end
abstract type AbstractModelNode <: AbstractStructureNode end
abstract type AbstractLinkingEdge  <: AbstractStructureEdge end
abstract type AbstractGraphSolver end

#Link Model
abstract type AbstractLinkConstraint <: JuMP.AbstractConstraint end
abstract type AbstractGraphConstraintRef end
abstract type AbstractLinkModel <: JuMP.AbstractModel end

include("linkmodel.jl")          #A JuMP extension model to manage GraphConstraints and LinkConstraints

include("modelnode.jl")          #ModelGraph nodes

include("modeledge.jl")          #ModelGraph edges

include("modelgraph.jl")         #The ModelGraph

include("macros.jl")             #@linkconstraint, @graphobjective

include("jumpgraph.jl")          #An aggregated JuMP model

include("solve.jl")              #Aggregate and solve with an MOI Solver
#
# include("solution.jl")         #SolutionGraph
#
# include("community_detection.jl")
#
include("partitioning/graph_projections.jl")  #Projections that facilitate graph analysis (partitioning and community detection)

include("partitioning/partition.jl")

include("partitioning/aggregation.jl")  #Aggregate pieces of a ModelGraph
#
#TODO: Update Metis to use JuliaSparse
# function __init__()
#     #@require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("extras/plots.jl")
#     @require Metis = "2679e427-3c69-5b7f-982b-ece356f1e94b" include("extras/metis.jl")
#     #@require CommunityDetection = "d427f087-d71a-5a1b-ace0-b93392eea9ff" include("extras/communities.jl")
# end


# if haskey(Pkg.installed(),"MPI")
# #External Solver Interfaces
#     include("solver_interfaces/wrapped_solvers.jl")
#
#     include("solver_interfaces/plasmoPipsNlpInterface.jl")
#     using .PlasmoPipsNlpInterface
# end

end
