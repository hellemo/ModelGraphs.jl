# 
# #IDEA: Go through variables and swap out references of NodeVariable to the underlying Variable
# function _splice_vrefs!(expr::Expr)  #var_map needs to map the node_model index to the new model variable
#     for i = 1:length(expr.args)
#         if typeof(expr.args[i]) == Expr
#             if expr.args[i].head != :ref             #keep calling _splice_nonlinear_variables! on the expression until it's a :ref. i.e. :(x[index])
#                 _splice_vrefs!(expr.args[i])
#             else  #it's a variable
#                 var_index = expr.args[i].args[2]     #this is the variable
#                 #var = :($(JuMP.VariableRef(model,var_index)))
#                 expr.args[i] = new_var               #replace :(x[index]) with a :(JuMP.Variable)
#             end
#         end
#     end
# end


#link variables can be added to a graph
macro linkvariable(graph,args...)
    code = quote
        @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
        JuMP.@variable($graph,($(args...)))
    end
    return esc(code)
end

macro link(graph,lvref,args...)
    code = quote
        for var in $(args...)
            link_variables!($graph,$lvref,var)
        end
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
        JuMP.@constraint($graph,($(args...)))   #this will call add_constraint(graph::ModelGraph)
    end
    return esc(code)
end

macro NLlinkconstraint(graph,args...)
    code = quote
        @assert isa($graph,AbstractModelGraph)  #Check the inputs are the correct types.  This needs to throw
        JuMP.@NLconstraint($graph,($(args...)))
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

#Create set of model nodes for multiple axes
#USAGE: @modelnodes(graph,t[1:24],x[1:10])
macro modelnodes(graph,args...)
end
