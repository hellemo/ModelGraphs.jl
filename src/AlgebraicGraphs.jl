module AlgebraicGraphs

using NestedHyperGraphs

using Requires
using LinearAlgebra
using DataStructures
using SparseArrays

using MathOptInterface
const MOI = MathOptInterface
const NHG = NestedHyperGraphs

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

#HyperPartition
HyperPartition,

#Solvers/Optimizers
AbstractGraphOptimizer,BendersOptimizer,DDOptimizer, ADMMOptimizer,

#Graph Functions
add_node!,getnode,getnodes,getlinkedges,getnumnodes,gethypergraph,

#Model functions
set_model,set_optimizer,reset_model,is_nodevariable,is_linked_variable,getmodel,has_model,
link_variables!,
#Aggregation
aggregate,aggregate!,

#solve handles
optimize!,bendersolve,dual_decomposition_solve,

#Solution management
nodevalue,nodedual,

#macros
@linkconstraint, @NLlinkconstraint, @linkvariable,

@NLnodeconstraint,

@masterconstraint,@NLmasterconstraint,

@graphobjective, @NLgraphobjective


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

include("aggregation.jl")          #An aggregated JuMP model

include("solve.jl")              #Aggregate and solve with an MOI Solver

# include("hyperpartition")

function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plots.jl")
end

end
