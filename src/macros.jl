macro linkvariable(graph,args...)
    code = quote
        @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
        JuMP.@variable($graph,($(args...)))
    end
    return esc(code)
end

#Constraints on link variables
macro masterconstraint(graph,args...)
    code = quote
        @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
        JuMP.@constraint($graph,($(args...)))    #this will call add_constraint(graph::ModelGraph)
    end
    return esc(code)
end

macro NLmasterconstraint(graph,args...)
    code = quote
        @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
        JuMP.@NLconstraint($(graph.mastermodel),($(args...)))  #link model extends @constraint macro
    end
    return esc(code)
end

#Node Constraints
#Wrap NLconstraint because NLconstraint extensions don't really work yet.  Easy to deprecate later.
macro NLnodeconstraint(node,args...)
    code = quote
        @assert isa($node,ModelNode)  #Check the inputs are the correct types.  This needs to throw
        JuMP.@NLconstraint((getmodel($node)),($(args...)))  #link model extends @constraint macro
    end
    return esc(code)
end


macro linkconstraint(graph,args...)
    code = quote
        @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
        JuMP.@constraint($graph,($(args...)))    #this will call add_constraint(graph::ModelGraph)
    end
    return esc(code)
end

macro NLlinkconstraint(graph,args...)
    code = quote
        @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
        JuMP.@NLconstraint($graph,($(args...)))  #link model extends @constraint macro
    end
    return esc(code)
end


#Graph Objectives
macro graphobjective(graph,args...)
    code = quote
        @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
        JuMP.@objective($graph,($(args...)))  #link model extends @constraint macro
    end
    return esc(code)
end

macro NLgraphobjective(graph,args...)
    code = quote
        @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
        JuMP.@NLobjective($graph,($(args...)))  #link model extends @constraint macro
    end
    return esc(code)
end
