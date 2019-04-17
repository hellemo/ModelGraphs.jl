"""
    @linkconstraint(graph,args...)
    macro for defining linkconstraints between nodes Link constraints are associated with nodes within their respective graph.
"""
macro linkconstraint(graph,args...)
        code = quote
            @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
            link_model = AlgebraicGraphs.getlinkmodel($graph)

            link_constraints = JuMP.@constraint(link_model,($(args...)))  #link model extends @constraint macro

            #TODO Check hypergraph implementation.  Fixed issues with slowness.
            AlgebraicGraphs.addlinkedges!($graph,link_constraints) #Go through each link constraint and add the appropriate edge
        end
        return esc(code)
end
