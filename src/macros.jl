"""
    @linkconstraint(graph,args...)
    macro for defining linkconstraints between nodes.
"""
macro linkconstraint(graph,args...)
        code = quote
            @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
            link_model = AlgebraicGraphs.getlinkmodel($graph)

            link_constraints = JuMP.@constraint(link_model,($(args...)))  #link model extends @constraint macro

            #TODO  Check the hypergraph implementation. I fixed the issues with slowness, but haven't tested it enough.
            #TODO: Create an add_constraint function for a ModelGraph.  It should add to the LinkModel and then add edges to the graph
            #AlgebraicGraphs.addlinkedges!($graph,link_constraints) #Go through each link constraint and add the appropriate edge
        end
        return esc(code)
end

#TODO: Finish macros
macro graphvariable
end

# A constraint between graph variables or graph variables and node variables
# NOTE: Hierarchy should be enforced.
macro graphconstraint
end

macro graphobjective
end
