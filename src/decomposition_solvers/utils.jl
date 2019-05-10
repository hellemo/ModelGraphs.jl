#set objective values to minimize for each node
function normalizegraph(graph::AbstractModelGraph)
    n = 1
    for node in getnodes(graph)
        m = getmodel(node)
        #Maximize --> Minimize the negative
        if JuMP.objective_sense == MOI.MAX_SENSE
            JuMP.set_objective_sense(m,MOI.MIN_SENSE)
            JuMP.set_objective_function(m,-1*JuMP.objective_function)
            n = -1
        end
    end
    setattribute(graph, :normalized, n)
    return n
end

function fix(var::JuMP.Variable,value::Real)
    setlowerbound(var,value)
    setupperbound(var,value)
end

#NOTE: I'm getting rid of the ModelTree object.  Everything is a graph now.  children nodes are nodes in a subgraph
"""
Checks if n1 is a child node of n2
"""
ischildnode(tree::ModelTree, n1::ModelNode, n2::ModelNode) = in(n2,in_neighbors(tree,n1))

function savenodeobjective(mf::JuMP.Model)
    g = mf.ext[:Graph]
    numnodes = length(getnodes(g))
    nodeindex = Dict("node$i" => i for i in 1:numnodes)
    nov = mf.ext[:nodeobj] = [AffExpr(0.0) for i in 1:numnodes]
    obj = mf.obj.aff
    for (i,var) in enumerate(obj.vars)
        coeff = obj.coeffs[i]
        varname = mf.colNames[var.col]
        nodename = varname[1:search(varname,'.')-1]
        index = nodeindex[nodename]
        push!(nov[index],coeff,var)
    end
end
