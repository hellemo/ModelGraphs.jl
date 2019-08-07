

macro linkvariable(graph,args...)
    code = quote
        @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
        @variable(graph,args...)
    end
    return esc(code)
end

macro masterconstraint(graph,args...)
    code = quote
        @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
        JuMP.@constraint(graph.master,($(args...)))    #this will call add_constraint(graph::ModelGraph)
    end
    return esc(code)
end

macro NLmasterconstraint(graph,args...)
    code = quote
        @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
        JuMP.@NLconstraint(graph.mastermodel,($(args...)))  #link model extends @constraint macro
    end
    return esc(code)
end

#Wrap NLconstraint because NLconstraint extensions don't work yet.  Easy to deprecate later.
macro NLnodeconstraint(node,args...)
end

# A constraint between graph variables or graph variables and node variables
"""
    @linkconstraint(graph,args...)
    macro for defining linkconstraints between nodes.
"""
macro linkconstraint(graph,args...)
        code = quote
            @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
            JuMP.@constraint(graph,($(args...)))    #this will call add_constraint(graph::ModelGraph)
        end
        return esc(code)
end

macro NLlinkconstraint(graph,args...)
        code = quote
            @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
            link_model = AlgebraicGraphs.getlinkmodel($graph)
            JuMP.@NLconstraint(link_model,($(args...)))  #link model extends @constraint macro
            #TODO  Check the hypergraph implementation. I fixed the issues with slowness, but haven't tested it enough.
        end
        return esc(code)
end

macro graphobjective
end

macro NLgraphobjective
end
