module AlgebraicGraphs

using NestedHyperGraphs

using Requires
using Distributed
using LinearAlgebra
using DataStructures
using SparseArrays

using MathOptInterface
const MOI = MathOptInterface

using JuMP
macro exportall(pkg)
    Expr(:export, names(JuMP)...)
end
@exportall JuMP

import JuMP: AbstractModel, AbstractConstraint, AbstractJuMPScalar, ConstraintRef
import Base: ==,show,print,string,getindex,copy

#Model Graph
export AbstractModelGraph, ModelGraph,

#Nodes and Edges,
ModelNode,LinkEdge,

#LinkVariable and LinkConstraint
LinkVariable, LinkConstraint, LinkVariableRef, LinkConstraintRef,

#ModelPartition
ModelPartition,

#Solvers
AbstractGraphOptimizer,BendersOptimizer,DualDecompositionOptimizer,

#Graph Functions
add_node!,getnode,getnodes,getlinkedges,getnumnodes,gethypergraph,

#Model functions
set_model,set_optimizer,reset_model,is_nodevariable,get_model,has_model,

#Aggregation
aggregate,aggregate!,

#solve handles
optimize!,bendersolve,dual_decomposition_solve,

#Solution management
nodevalue,nodedual,

#macros
@linkconstraint,@linkvariable,@NLlinkconstraint,@graphobjective

#Abstract Types
abstract type AbstractModelGraph <: JuMP.AbstractModel end
abstract type AbstractLinkEdge end
abstract type AbstractLinkConstraintRef end
abstract type AbstractLinkVariableRef <: JuMP.AbstractVariableRef end
abstract type AbstractGraphOptimizer end
abstract type AbstractLinkConstraint <: JuMP.AbstractConstraint end

include("modelnode.jl")          #ModelGraph nodes

include("linkedge.jl")          #ModelGraph edges

include("modelgraph.jl")         #The ModelGraph

#include("nlp_extension.jl")

# include("macros.jl")             #@linkconstraint, @graphobjective

# include("modelpartition")
#
# include("aggregate.jl")          #An aggregated JuMP model
#
# include("solve.jl")              #Aggregate and solve with an MOI Solver
#
# include("solutiongraph.jl")         #SolutionGraph
#
# include("community_detection.jl")
#
# include("partitioning/graph_projections.jl")  #Projections that facilitate graph analysis (partitioning and community detection)
#
# include("package_extensions/metis.jl")
#
# include("partitioning/partition.jl")
#
# include("partitioning/aggregation.jl")  #Aggregate pieces of a ModelGraph


# #Decomposition-based solvers
# include("decomposition_solvers/utils.jl")
#
# include("decomposition_solvers/solution.jl")
#
# include("decomposition_solvers/lagrange/dual_decomposition.jl")


# function __init__()
#     #@require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("extras/plots.jl")
#     @require Metis = "2679e427-3c69-5b7f-982b-ece356f1e94b" include("extras/metis.jl")
#     #@require CommunityDetection = "d427f087-d71a-5a1b-ace0-b93392eea9ff" include("extras/communities.jl")
# end

end
