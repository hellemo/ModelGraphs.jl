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

#Solvers/Optimizers
AbstractGraphOptimizer,BendersOptimizer,DualDecompositionOptimizer,

#Graph Functions
add_node!,getnode,getnodes,getlinkedges,getnumnodes,gethypergraph,

#Model functions
set_model,set_optimizer,reset_model,is_nodevariable,getmodel,has_model,

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

include("modelnode.jl")

include("linkedge.jl")

include("modelgraph.jl")

include("nlp_extension.jl")

include("macros.jl")

# include("hyperpartition")

# include("aggregation.jl")          #An aggregated JuMP model

# include("solve.jl")              #Aggregate and solve with an MOI Solver
#
# include("solution.jl")         #SolutionGraph


# #Decomposition-based solvers
# include("decomposition_solvers/utils.jl")
#
# include("decomposition_solvers/solution.jl")
#
# include("decomposition_solvers/lagrange/dual_decomposition.jl")
# include("decomposition_solvers/benders/dual_decomposition.jl")


# function __init__()
#     #@require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("extras/plots.jl")
# end

end
