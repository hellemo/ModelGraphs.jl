
function get_link_eq_constraints(graph::ModelGraph)
    link_cons = LinkConstraint[]
    for link_con in getlinkconstraints(graph)
        if isa(link_con.set,MOI.EqualTo)
            push!(link_cons,link_con)
        end
    end
    return link_cons
end

function get_link_ineq_constraints(graph::ModelGraph)
    link_cons = LinkConstraint[]
    for link_con in getlinkconstraints(graph)
        if !(isa(link_con.set,MOI.EqualTo))
            push!(link_cons,link_con)  #It must be an inequality constraint
        end
    end
    return link_cons
end

function prepare_link_matrix(link_constraints::Vector{LinkConstraint})
    #link_constraints = graph.linkmodel.linkconstraints #index => constraint
    link_map = Dict()
    link_vars = []
    n_link_vars = 0
    I = Int64[]
    J = Int64[]
    V = Float64[]
    b = Float64[]
    i = 0
    for link in link_constraints  #Just enumerate here
        i += 1
        _unpack_constraint!(link,i,I,J,V,b,link_vars,link_map)
    end
    A = sparse(I,J,V)
    return A,b,link_vars,link_map
end

function _unpack_constraint!(constraint::LinkConstraint{GenericAffExpr{Float64,VariableRef},MOI.EqualTo{Float64}},row::Integer,I::Vector,J::Vector,V::Vector,b::Vector,link_vars::Vector,link_map::Dict)
    push!(b,constraint.set.value)
    for (var,coeff) in constraint.func.terms
        #find a new link variable
        if !(var in link_vars)
            j = length(link_vars) + 1
            link_map[var] = j   #link index of variable var
            push!(link_vars,var)
        else
            #get correct variable index
            j = link_map[var]
        end
        push!(I,row)
        push!(J,j)
        push!(V,coeff)
    end
end

#NOTE: Might make this into a macro
function _unpack_constraint!(constraint::LinkConstraint{GenericAffExpr{Float64,VariableRef},MOI.LessThan{Float64}},row::Integer,I::Vector,J::Vector,V::Vector,b::Vector,link_vars::Vector,link_map::Dict)
    push!(b,constraint.set.upper)
    for (var,coeff) in constraint.func.terms
        #find a new link variable
        if !(var in link_vars)
            j = length(link_vars) + 1
            link_map[var] = j   #link index of variable var
            push!(link_vars,var)
        else
            #get correct variable index
            j = link_map[var]
        end
        push!(I,row)
        push!(J,j)
        push!(V,coeff)
    end
end

function _unpack_constraint!(constraint::LinkConstraint{GenericAffExpr{Float64,VariableRef},MOI.GreaterThan{Float64}},row::Integer,I::Vector,J::Vector,V::Vector,b::Vector,link_vars::Vector,link_map::Dict)
    push!(b,-1*constraint.set.lower)
    for (var,coeff) in constraint.func.terms
        #find a new link variable
        if !(var in link_vars)
            j = length(link_vars) + 1
            link_map[var] = j   #link index of variable var
            push!(link_vars,var)
        else
            #get correct variable index
            j = link_map[var]
        end
        push!(I,row)
        push!(J,j)
        push!(V,-1*coeff)
    end
end
#TODO: Interval

# function prepare_link_ineq_matrix(graph::ModelGraph)
#     #link_constraints = graph.linkmodel.linkconstraints #index => constraint
#     link_map = Dict()
#     link_vars = []
#     n_link_vars = 0
#     I = Int64[]
#     J = Int64[]
#     V = Float64[]
#     b = Float64[]
#     i = 0
#     for link_con in getlinkconstraints(graph)
#         #If it's an equality constraint
#
#         if isa(link_con.set,MOI.LessThan)  #Good to go
#             i += 1
#
#             push!(b,link_con.set.upper)
#             #Don't need to modify coefficients
#             for (var,coeff) in link_con.func.terms
#                 #find a new link variable
#                 if !(var in keys(link_vars))
#                     n_link_vars += 1
#                     j = n_link_vars
#                     link_map[var] = j   #link index of variable var
#                     push!(link_vars,var)
#                 else
#                     #get correct variable index
#                     j = link_vars[var]
#                 end
#
#                 push!(I,i)
#                 push!(J,j)
#                 push!(V,coeff)
#             end
#
#         elseif isa(link_con.set,MOI.GreaterThan)
#             i += 1
#             push!(b,-1*link_con.set.lower)
#             #Don't need to modify coefficients
#             for (var,coeff) in link_con.func.terms
#                 #find a new link variable
#                 if !(var in keys(link_vars))
#                     n_link_vars += 1
#                     j = n_link_vars
#                     link_map[var] = j   #link index of variable var
#                     push!(link_vars,var)
#                 else
#                     #get correct variable index
#                     j = link_vars[var]
#                 end
#
#                 push!(I,i)
#                 push!(J,j)
#                 push!(V,-1*coeff)
#             end
#
#         elseif isa(link_con.set,MOI.Interval)
#             error("Lagrange Solver does not yet support Interval constraints")
#         end
#     end
#
#     Pi = sparse(I,J,V)
#     return Pi,link_vars,b
# end



# #set objective values to minimize for each node
# function normalizegraph(graph::AbstractModelGraph)
#     n = 1
#     for node in getnodes(graph)
#         m = getmodel(node)
#         #Maximize --> Minimize the negative
#         if JuMP.objective_sense == MOI.MAX_SENSE
#             JuMP.set_objective_sense(m,MOI.MIN_SENSE)
#             JuMP.set_objective_function(m,-1*JuMP.objective_function)
#             n = -1
#         end
#     end
#     setattribute(graph, :normalized, n)
#     return n
# end
#
# function fix(var::JuMP.Variable,value::Real)
#     setlowerbound(var,value)
#     setupperbound(var,value)
# end
#
# #NOTE: I'm getting rid of the ModelTree object.  Everything is a graph now.  children nodes are nodes in a subgraph
# """
# Checks if n1 is a child node of n2
# """
# ischildnode(tree::ModelTree, n1::ModelNode, n2::ModelNode) = in(n2,in_neighbors(tree,n1))
#
# function savenodeobjective(mf::JuMP.Model)
#     g = mf.ext[:Graph]
#     numnodes = length(getnodes(g))
#     nodeindex = Dict("node$i" => i for i in 1:numnodes)
#     nov = mf.ext[:nodeobj] = [AffExpr(0.0) for i in 1:numnodes]
#     obj = mf.obj.aff
#     for (i,var) in enumerate(obj.vars)
#         coeff = obj.coeffs[i]
#         varname = mf.colNames[var.col]
#         nodename = varname[1:search(varname,'.')-1]
#         index = nodeindex[nodename]
#         push!(nov[index],coeff,var)
#     end
# end
