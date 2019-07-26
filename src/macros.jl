"""
    @linkconstraint(graph,args...)
    macro for defining linkconstraints between nodes.
"""
macro linkconstraint(graph,args...)
        code = quote
            @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
            link_model = AlgebraicGraphs.getlinkmodel($graph)
            JuMP.@constraint(link_model,($(args...)))  #link model extends @constraint macro
            #TODO  Check the hypergraph implementation. I fixed the issues with slowness, but haven't tested it enough.
        end
        return esc(code)
end

#TODO: Finish macros
macro linkvariable
end

#Wrap NLconstraint because NLconstraint extensions don't work yet.  Easy to deprecate later.
macro NLnodeconstraint(node,args...)
end
# A constraint between graph variables or graph variables and node variables
# NOTE: Hierarchy should be enforced.
macro graphconstraint
end

macro graphobjective
end
