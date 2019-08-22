ModelGraphs.jl is currently a work in progress.  The first release should happen soon.

# ModelGraphs.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jalving.github.io/AlgebraicGraphs.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jalving.github.io/AlgebraicGraphs.jl/dev)
[![Build Status](https://travis-ci.com/jalving/AlgebraicGraphs.jl.svg?branch=master)](https://travis-ci.com/jalving/AlgebraicGraphs.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jalving/AlgebraicGraphs.jl?svg=true)](https://ci.appveyor.com/project/jalving/AlgebraicGraphs-jl)
[![Codecov](https://codecov.io/gh/jalving/AlgebraicGraphs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jalving/AlgebraicGraphs.jl)
[![Coveralls](https://coveralls.io/repos/github/jalving/AlgebraicGraphs.jl/badge.svg?branch=master)](https://coveralls.io/github/jalving/AlgebraicGraphs.jl?branch=master)

ModelGraphs.jl is a package that provides a graph-based approach to algebraic modeling.  It makes it possible to model optimization structures in a modular way wherein component models (model nodes) can be connected
in a graph using link constraints and link variables.  This convenient modeling paradigm also facilitates modeling problems in a structured way to use distributed optimization solvers that can exploit block diagonal forms.  A model graph also maintains an underlying hypergraph representation which can be used through a generic partitioning interface to aggregate complex optimization structures into standard structured forms.

## Installation

ModelGraphs.jl currently requires an unregistered dependency NestedHyperGraphs.jl.  You will need to install the dependency and then ModelGraphs.jl as follows:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/jalving/NestedHyperGraphs.jl.git"))
Pkg.add(PackageSpec(url="https://github.com/jalving/ModelGraphs.jl.git"))
```


## Simple Example

```julia
using ModelGraphs
using Ipopt

graph = ModelGraph()
optimizer = with_optimizer(Ipopt.Optimizer)

#Add nodes to a ModelGraph
n1 = add_node!(graph)
n2 = add_node!(graph)

@variable(n1,0 <= x <= 2)
@variable(n1,0 <= y <= 3)
@constraint(n1,x+y >= 4)

@variable(n2,x)
@NLnodeconstraint(n2,ref,exp(x) >= 2)

#Link constraints take the same expressions as the JuMP @constraint macro
@linkconstraint(graph,n1[:x] == n2[:x])
@graphobjective(graph,Min,n1[:y] + n2[:x])

optimize!(graph,optimizer)

println("n1[:x]= ",nodevalue(n1[:x]))
println("n1[:y]= ",nodevalue(n1[:y]))
println("n2[:x]= ",nodevalue(n2[:x]))
println("objective = ", objective_value(graph))
```
