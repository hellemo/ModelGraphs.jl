# ModelGraphs

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jalving.github.io/AlgebraicGraphs.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jalving.github.io/AlgebraicGraphs.jl/dev)
[![Build Status](https://travis-ci.com/jalving/AlgebraicGraphs.jl.svg?branch=master)](https://travis-ci.com/jalving/AlgebraicGraphs.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jalving/AlgebraicGraphs.jl?svg=true)](https://ci.appveyor.com/project/jalving/AlgebraicGraphs-jl)
[![Codecov](https://codecov.io/gh/jalving/AlgebraicGraphs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jalving/AlgebraicGraphs.jl)
[![Coveralls](https://coveralls.io/repos/github/jalving/AlgebraicGraphs.jl/badge.svg?branch=master)](https://coveralls.io/github/jalving/AlgebraicGraphs.jl?branch=master)

ModelGraphs.jl is a package that provides a graph-based approach to algebraic modeling.  It makes it possible to model optimization structures in a modular way wherein component models (model nodes) can be connected
in a graph using link constraints and link variables.  This convenient modeling paradigm also facilitates modeling problems in a structured way to use distributed optimization solvers that can exploit block diagonal forms.  A model graph also maintains an underlying hypergraph representation which can be used through a generic partitioning interface to aggregate complex optimization structures into standard structured forms.
