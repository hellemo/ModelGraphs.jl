#####################################################
#  Link Constraint
#  A linear constraint between JuMP Models (nodes).  Link constraints can be equality or inequality.
#####################################################
struct LinkConstraint{F <: AbstractJuMPScalar,S <: MOI.AbstractScalarSet} <: AbstractLinkConstraint
    func::F
    set::S
    node_indices::Vector{Int}  #Custom link constraint data.  Indices of nodes this LinkConstraint connects.
end

getnode(m::JuMP.Model) = m.ext[:node]

function getnode(var::JuMP.VariableRef)
    if haskey(var.model.ext,:node)
        return getnode(var.model)
    else
        error("variable $var does not belong to a node.  If you're trying to create a linkconstraint, make sure
        the owning model has been set to a node.")
    end
end

function LinkConstraint(con::JuMP.ScalarConstraint)
    node_indices = [getnode(var) for var in keys(con.func.terms)]
    return LinkConstraint(con.func,con.set,node_indices)
end

function JuMP.add_constraint(m::LinkModel, con::JuMP.ScalarConstraint, name::String="")
    m.graph_constraint_index += 1
    cref = GraphConstraintRef(m, m.graph_constraint_index)
    link_con = LinkConstraint(con)
    m.linkconstraints[cref.idx] = link_con
    JuMP.set_name(cref, name)
    return cref
end

function JuMP.add_constraint(m::LinkModel, con::JuMP.AbstractConstraint, name::String="")
    error("Link Models only support graph constraints and link constraints")
end

is_simplelinkconstr(con::LinkConstraint) = getnumnodes(con) == 2 ? true : false
is_hyperlinkconstr(con::LinkConstraint) = getnumnodes(con) > 2 ? true : false

jump_function(constraint::LinkConstraint) = constraint.func
moi_set(constraint::LinkConstraint) = constraint.set
shape(::LinkConstraint) = JuMP.ScalarShape()

"Add a single link-constraint to the ModelGraph"
function addlinkconstraint(graph::AbstractModelGraph,linkcon::AbstractLinkConstraint)
    isa(con,JuMP.LinearConstraint) || throw(error("Link constraints must be linear.  If you're trying to add quadtratic or nonlinear links, try creating duplicate variables and linking those"))
    ref = JuMP.addconstraint(graph.linkmodel,con)
    link_edge = add_edge!(graph,ref)  #adds edge and a contraint reference to all objects involved in the constraint
    return link_edge
end

#NOTE Figure out a good way to use containers here instead of making arrays
"Add a vector of link-constraints to the ModelGraph"
function addlinkconstraint(graph::AbstractModelGraph,linkcons::Array{AbstractConstraint,T}) where T
    #NOTE I don't know why I wrote these two lines anymore
    #array_type = typeof(linkcons)   #get the array type
    #array_type.parameters.length > 1 ? linkcons = vec(linkcons) : nothing   #flatten out the constraints into a single vector
    linkcons = vec(linkcons)

    #Check all of the constraints before I add one to the graph
    for con in linkcons
        vars = con.terms.vars
        nodes = unique([getnode(var) for var in vars])
        all(node->node in getnodes(graph),nodes) ? nothing : error("the linkconstraint: $con contains variables that don't belong to the graph: $graph")
    end

    #Now add the constraints
    for con in linkcons
        addlinkconstraint(graph,con)
    end
end

"""
getlinkconstraints(graph::AbstractModelGraph)

Return Array of all LinkConstraints in the ModelGraph graph
"""
getlinkconstraints(model::AbstractModelGraph) = getlinkconstraints(model.linkmodel)

"""
getsimplelinkconstraints(model::AbstractModelGraph)

Retrieve link-constraints that only connect two nodes"
"""
getsimplelinkconstraints(model::AbstractModelGraph) = getsimplelinkconstraints(model.linkmodel)


"""
gethyperlinkconstraints(model::AbstractModelGraph)

Retrieve link-constraints that connect three or more nodes"
"""
gethyperlinkconstraints(model::AbstractModelGraph) = gethyperlinkconstraints(model.linkmodel)

"""
get_all_linkconstraints(graph::AbstractModelGraph)

Get a list containing every link constraint in the graph, including its subgraphs
"""
function get_all_linkconstraints(graph::AbstractModelGraph)
    links = []
    for subgraph in getsubgraphlist(graph)
        append!(links,getlinkconstraints(subgraph))
    end
    append!(links,getlinkconstraints(graph))
    return links
end
