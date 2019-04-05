#IDEA : If solving with JuMP: Create the corresponding JuMP model and return the solution to the graph

# struct ReferenceMap
#     model::Model
#     index_map::MOIU.IndexMap
# end
# function Base.getindex(reference_map::ReferenceMap, vref::VariableRef)
#     return JuMP.VariableRef(reference_map.model,
#                        reference_map.index_map[index(vref)])
# end
# function Base.getindex(reference_map::ReferenceMap, cref::ConstraintRef)
#     return JuMP.ConstraintRef(reference_map.model,
#                          reference_map.index_map[index(cref)],
#                          cref.shape)
# end
# Base.broadcastable(reference_map::ReferenceMap) = Ref(reference_map)


mutable struct JuMPGraph <: AbstractModelGraph
    hypergraph::StructuredHyperGraph
    graphvariables::Vector{JuMP.VariableRef}
    graphconstraints::Vector{JuMP.ConstraintRef}
    linkconstraints::Vector{JuMP.ConstraintRef}
    #referencemap::JuMP.ReferenceMap  #Map subproblems variable and constraint references
end
JuMPGraph() = JuMPGraph(BasePlasmoGraph(HyperGraph),ConstraintRef[])
JuMPGraph(basegraph::BasePlasmoGraph) = JuMPGraph(basegraph,ConstraintRef[])

mutable struct JuMPNode <: AbstractModelNode
    node::StructureNode
    obj_dict::Dict{Symbol,Any}
    variablelist::Vector{VariableRef}
    constraintlist::Vector{ConstraintRef}
    #index_map::MOIU.IndexMap                        #Map of variable and constraint indices from JuMP node to aggregated model
    #linear index of node variable in new jump graph model to the original index of the node model
end
create_node(graph::JuMPGraph) = JuMPNode(StructureNode(),0,Dict{Symbol,AbstractJuMPScalar}(),AbstractJuMPScalar[],ConstraintRef[],Dict{Int,Int}())
#hasmodel(node::JuMPNode) = throw(error("JuMP nodes are simple references to original ModelNodes.  Did you mean to check a ModelNode?"))

#Has constraint references for link constraints
mutable struct JuMPEdge <: AbstractLinkingEdge
    edge::StructureEdge
    linkconstraints::Vector{ConstraintRef}  #indices in JuMP model of linkconstraints for this edge
end
create_edge(graph::JuMPGraph) = JuMPEdge(BasePlasmoEdge(),ConstraintRef[])

#Construct a structured model, but roll it all into one JuMP model (this is how we solve with JuMP accessible solvers)
function JuMPGraphModel()
    m = JuMP.Model()
    m.ext[:Graph] = JuMPGraph()
    return m
end
is_graphmodel(m::JuMP.Model) = haskey(m.ext,:Graph) ? true : false  #check if the model is a graph model

# Should be defined by base type
# #Add nodes and edges to graph models.  These are used for model instantiation from a graph
function StructureGraphs.add_node!(m::JuMP.Model; index = nv(getgraph(m).hypergraph)+1)
    is_graphmodel(m) || error("Can only add nodes to graph models")
    node = create_node(getgraph(m))
    add_node!(getgraph(m),node,index = index)
    return node
end

function StructureGraphs.add_edge!(m::Model,nodes::JuMPNode...)
    is_graphmodel(m) || error("Can only add edges to graph models")
    edge = add_edge!(getgraph(m),nodes...)
    return edge
end

#Define all of the JuMP model extension functions
getgraph(m::Model) = haskey(m.ext, :Graph) ? m.ext[:Graph] : error("Model is not a graph model")
StructureGraphs.getnodes(m::Model) = getnodes(getgraph(m))
StructureGraphs.getedges(m::Model) = getedges(getgraph(m))

StructureGraphs.getnode(m::Model,id::Integer) = getnode(getgraph(m),id)  #Grab from the highest level graph if not specified
StructureGraphs.getnode(m::Model,sid::Integer,nid::Integer) = getnode(getgraph(m).subgraphlist[sid],nid)

StructureGraphs.getedge(m::Model,id::LightGraphs.AbstractEdge) = getedge(getgraph(m))[id]
StructureGraphs.getedge(m::Model,sid::Integer,eid::LightGraphs.AbstractEdge) = getedge(getgraph(m).subgraphlist[sid],eid)

JuMP.getobjective(node::JuMPNode) = node.objective
JuMP.getobjectivevalue(node::JuMPNode) = getvalue(node.objective)

getnodevariablemap(node::JuMPNode) = node.variablemap
getnodevariables(node::JuMPNode) = node.variablelist
getnodevariable(node::JuMPNode,index::Integer) = node.variablelist[index]
getnodeconstraints(node::JuMPNode) = node.constraintlist
JuMP.num_variables(node::JuMPNode) = length(node.variablelist)
#get node variables
getindex(node::JuMPNode,s::Symbol) = node.variablemap[s]

#get all of the link constraints from a JuMP model
#Retrieves constraint indices that are link constraints
function getlinkconstraints(m::JuMP.Model)
    is_graphmodel(m) || error("link constraints are only available on graph models")
    return getgraph(m).linkconstraints
end


#reference_maps = Dict{JuMPNode,JuMP.ReferenceMap}()
#var_maps = Dict{JuMPNode,Dict}()
#index_maps = Dict{JuMPNode,IndexMap}()

#Create a single JuMP model from a ModelGraph
function create_jump_graph_model(model_graph::AbstractModelGraph;add_node_objectives = hasobjective(model_graph))
    jump_graph_model = JuMPGraphModel()

    jump_graph = copy_graph_to(model_graph,JuMPGraph)  #Create empty JuMP graph with nodes and edges
    jump_graph_model.ext[:Graph] = jump_graph

    # if add_node_objectives
    #     graph_objective = sum(objective(node) for node )
    ref_map = JuMP.ReferenceMap(jump_graph_model,MOIU.IndexMap())
    #COPY NODE MODELS INTO AGGREGATED MODEL
    for model_node in getnodes(model_graph)             #for each node in the model graph
        nodeindex = getindex(model_graph,model_node)
        jump_node = getnode(jump_graph,nodeindex)
        node_ref_map = _buildnodemodel!(jump_graph_model,jump_node,model_node)
        reference_maps[jump_node] = node_ref_map

        #Add node objective
        if add_node_objectives
            node_objective = copy_objective(model_node)
            add_to_objective!(jump_graph_model)
        end
    end

    #linkconstraint = LinkConstraint(link)
    #indexmap = Dict() #{node variable => new jump model variable index} Need index of node variables to flat model variables
    #vars = keys(linkconstraint.func.terms)
    #LINK CONSTRAINTS
    #inspect the link constraints, and map them to variables within flat model
    for linkconstraint in get_all_linkconstraints(model_graph)
        new_constraint = copy_constraint(linkconstraint,reference_map)
        JuMP.add_constraint(m,new_constraint)
    end

    #GRAPH VARIABLES
    for graph_variable in model_graph.graphvariables
    end

    #GRAPH CONSTRAINTS
    for graph_constraint in model_graph.graphconstraints
    end

    #GRAPH OBJECTIVE
    #sum the objectives by default
    has_nonlinear_obj = false   #check if any nodes have nonlinear objectives
    for node in getnodes(model_graph)
        node_model = getmodel(node)
        nlp = node_model.nlpdata
        if nlp != nothing && nlp.nlobj != nothing
            has_nonlinear_obj = true
            break
        end
    end

    if has_nonlinear_obj  == false      #just sum linear or quadtratic objectives
        @objective(jump_model,Min,sum(getobjective(node) for node in getnodes(jump_graph)))

    elseif has_nonlinear_obj == true  #build up the objective expression and splice in variables.  Cast all objectives as nonlinear
        obj = :(0)
        #for (id,node) in getnodesandedges(flat_graph)
        for node in getnodes(jump_graph)
            id = getindex(jump_graph,node)
            node_model = getmodel(getnode(model_graph,id))
            getobjectivesense(node_model) == :Min ? sense = 1 : sense = -1

            #nlp = node_model.nlpdata
            #if nlp == nothing# || (nlp != nothing && nlp.nlobj == nothing) #cast the problem as nonlinear
            d = JuMP.NLPEvaluator(node_model)
            MOI.initialize(d,[:ExprGraph])
            node_obj = MOI.objective_expr(d)
            #_splice_nonlinear_variables!(node_obj,var_maps[node])
            _splice_nonlinear_variables!(node_obj,ref_map)

            node_obj = Expr(:call,:*,:($sense),node_obj)
            obj = Expr(:call,:+,obj,node_obj)
        end
        #println(obj)
        JuMP.set_NL_objective(jump_model, :Min, obj)
    end
    return jump_model
end



# elseif  nlp != nothing# && nlp.nlobj == nothing
#     d = JuMP.NLPEvaluator(node_model)
#     MathProgBase.initialize(d,[:ExprGraph])
#     node_obj = MathProgBase.obj_expr(d)
#     _splice_nonlinear_variables!(node_obj,var_maps[node])
# end
#node_obj = getnodeobjective(node)

#Useful utility for copying object_data from node models to a the aggregated model
# struct ReferenceMap
#     model::Model
#     index_map::Dict{Int,Int}
# end
# function Base.getindex(reference_map::ReferenceMap, vref::VariableRef)
#     return JuMP.VariableRef(reference_map.model,reference_map.index_map[index(vref)])
# end
# function Base.getindex(reference_map::ReferenceMap, cref::ConstraintRef)
#     return JuMP.ConstraintRef(reference_map.model,reference_map.index_map[index(cref)],cref.shape)
# end
# Base.broadcastable(reference_map::ReferenceMap) = Ref(reference_map)
#If I end up trying to this the JuMP way.

function copy_constraint_func(func::JuMP.GenericAffExpr{Float64,V <: AbstractVariableRef},ref_map::ReferenceMap)
    terms = func.terms
    new_terms = OrderedDict([(ref_map[var_ref],coeff) for (var_ref,coeff) in terms])
    new_func = JuMP.GenericAffExpr{Float64,typeof(V)}()
    new_func.terms = new_terms
    new_func.constant = func.constant
    return new_func
end

function copy_constraint_func(func::JuMP.GenericQuadExpr,ref_map::ReferenceMap)
    new_aff = copy_constraint_func(func.aff)
    new_terms = OrderedDict([(ref_map[var_ref],coeff) for (var_ref,coeff) in terms])
    new_func = JuMP.GenericQuadExpr{Float64,typeof(V)}()
    new_func.terms = new_terms
    new_func.aff = new_aff
    new_func.constant = func.constant
    return new_func
end

function copy_constraint_func(func::JuMP.VariableRef,ref_map::ReferenceMap)
    new_func = ref_map[func]
    return new_func
end

function copy_constraint(constraint::JuMP.ScalarConstraint,ref_map::ReferenceMap)
    new_func = copy_constraint_func(constraint.func,ref_map)
    new_con = JuMP.ScalarConstraint(new_func,constraint.set)
    return new_con
end

#NOTE Figure out whether I can do a broadcast for the array
function copy_constraint(constraints::JuMP.VectorConstraint,ref_map::ReferenceMap)
    new_funcs = [copy_constraint_func(con.func) for con in constraints]
    new_con = JuMP.VectorConstraint(new_func,constraint.set,constraint.shape)
    return new_con
end

function copy_constraint(constraint::GraphScalarConstraint)
    new_func = copy_constraint_func(constraint.func)
end

#Function to build a node model for a flat graph model
function _buildnodemodel!(m::Model,jump_node::JuMPNode,model_node::ModelNode)
    node_model = getmodel(model_node)

    if mode(model) == DIRECT
        error("Cannot copy a node model in `DIRECT` mode. Use the `Model` ",
              "constructor instead of the `direct_model` constructor to be ",
              "able to aggregate into a new JuMP Model.")
    end

    #Copy model_node into aggregate model
    #NOTE This might be another way to do the copy using the MOI interface
    #index_map = copy_node_model(m,node_model,true)

    #Copy object data (Containers)
    #See JuMP copy function for how this works
    #reference_map = ReferenceMap()

    #COPY VARIABLES
    index_map = MOIU.IndexMap()
    for var in JuMP.all_variables(node_model)
        new_x = JuMP.@variable(m)  #create an anonymous variable
        i = var.index              #This is an MOI Index

        #variable_map[i] = new_x                                                 #map the index of the model_node variable to the new variable in the jump_node
        index_map[i] = new_x.index                                               #map of model node variable index to the aggregated model index
                                     #Update master model variable dictionary
        #push!(jump_node.variablelist,new_x)
        var_name = JuMP.name(var)
        new_name = "$(getlabel(jump_node))$(getindex(getgraph(m),jump_node))."*var_name

        JuMP.set_name(x,new_name)
        JuMP.set_value(new_x,JuMP.value(var))
    end
    reference_map = ReferenceMap(node_model,index_map)


    #COPY CONSTRAINTS
    #NOTE
    #Option 1
    #Use JuMP and check if I have ScalarConstraint or VectorConstraint and use my index map to create new constraints

    #Option 2:
    #Use MOI to copy which would hit all of these constraints.  See MOI.copy_to for ideas about how to do this.
    constraint_types = JuMP.list_of_constraint_types(node_model)
    for (func,set) in constraint_types
        constraint_refs = JuMP.all_constraints(model, func, set)
        for constraint_ref in constraint_refs
            constraint = JuMP.constraint_object(contraint_ref)
            new_constraint = copy_constraint(constraint,index_map)
            JuMP.add_constraint(new_model,new_constraint)
        end
    end

    #COPY OBJECT DATA (JUMP CONTAINERS)
    for (name, value) in object_dictionary(node_model)
        jump_node.obj_dict[name] = getindex.(reference_map, value)
    end

    #COPY NONLINEAR CONSTRAINTS
    if node_model.nlp !== nothing
        d = JuMP.NLPEvaluator(node_model)           #Get the NLP evaluator object.  Initialize the expression graph
        MOI.initialize(d,[:ExprGraph])
        for nlconstr in node_model.nlpdata.nlconstr

        for i = start_index:num_cons
            expr = MOI.constraint_expr(d,i)  #this returns a julia expression
            _splice_nonlinear_variables!(expr,ref_map)              #splice the variables from var_map into the expression
            con = JuMP.add_NL_constraint(m,expr)    #raw expression input for non-linear constraint
            push!(jump_node.constraintlist,con)
        end
    end

    #OBJECTIVE FUNCTION
    JuMP.objectivesense(node_model) == :Min ? sense = 1 : sense = -1

    #LINEAR OR QUADRATIC OBJECTIVE

    nlp = node_model.nlpdata
    if nlp == nothing  || (nlp !== nothing && nlp.nlobj == nothing)
        #Get the linear terms
        t = []
        for terms in linearterms(node_model.obj.aff)
            push!(t,terms)
        end
        #Get the quadratic terms
        qcoeffs = node_model.obj.qcoeffs
        qvars1 = node_model.obj.qvars1
        qvars2 = node_model.obj.qvars2
        obj = @objective(m,Min,sense*(sum(qcoeffs[i]*var_map[linearindex(qvars1[i])]*var_map[linearindex(qvars2[i])] for i = 1:length(qcoeffs)) +
        sum(t[i][1]*var_map[linearindex(t[i][2])] for i = 1:length(t)) + node_model.obj.aff.constant))
        jump_node.objective = m.obj

    #NONLINEAR OBJECTIVE
    elseif nlp != nothing && nlp.nlobj != nothing
        obj = MathProgBase.obj_expr(d)
        _splice_nonlinear_variables!(obj,var_map)
        obj = Expr(:call,:*,:($sense),obj)
        jump_node.objective = m.obj
    end
    return m,var_map
end

#splice variables into a constraint expression
function _splice_nonlinear_variables!(expr::Expr,var_map::Dict)
    for i = 1:length(expr.args)
        if typeof(expr.args[i]) == Expr
            if expr.args[i].head != :ref   #keep calling _splicevars! on the expression until it's a :ref. i.e. :(x[index])
                _splice_nonlinear_variables!(expr.args[i],var_map)
            else  #it's a variable
                var_index = expr.args[i].args[2]     #this is the actual index in x[1], x[2], etc...
                new_var = :($(var_map[var_index]))   #get the JuMP variable from var_map using the index
                expr.args[i] = new_var               #replace :(x[index]) with a :(JuMP.Variable)
            end
        end
    end
end

#Create a JuMP model and solve with a MPB compliant solver
buildjumpmodel!(graph::AbstractModelGraph) = graph.serial_model = create_jump_graph_model(graph)
function jump_solve(graph::AbstractModelGraph;scale = 1.0,kwargs...)
    println("Aggregating Models...")
    m_flat = buildjumpmodel!(graph)
    println("Finished Creating JuMP Model")
    m_flat.solver = graph.linkmodel.solver
    m_flat.obj = scale*m_flat.obj
    status = JuMP.solve(m_flat;kwargs...)
    if status == :Optimal
        setsolution(getgraph(m_flat),graph)                       #Now get our solution data back into the original ModelGraph
        _setobjectivevalue(graph,JuMP.getobjectivevalue(m_flat))  #Set the graph objective value for easy access
    end
    return status
end

#check if graph has a mathprogbase solver
#TODO Remove scale argument.  Allow direct interface with the graph objective function
function JuMP.solve(graph::AbstractModelGraph;scale = 1.0,kwargs...)
    if isa(getsolver(graph),AbstractMathProgSolver)
        status = jump_solve(graph,scale = scale,kwargs...)
    elseif isa(getsolver(graph),AbstractPlasmoSolver)
        status = solve(graph,getsolver(graph))
    else
        throw(error("Given solver not recognized"))
    end
    return status
end



#copy the solution from one graph to another where nodes and variables match
function setsolution(graph1::AbstractModelGraph,graph2::AbstractModelGraph)
    for node in getnodes(graph1)
        index = getindex(graph1,node)
        node2 = getnode(graph2,index)       #get the corresponding node or edge in graph2
        for i = 1:num_var(node)
            node1_var = getnodevariable(node,i)
            node2_var = getnodevariable(node2,i)
            setvalue(node2_var,getvalue(node1_var))
        end
    end
    #TODO Set dual values for linear constraints
end

#IDEA: Use MOI to create a backend copy and then fill in the missing JuMP information
function my_copy_to(aggregate_model::JuMP.Model, src_model::JuMP.Model, copy_names::Bool)

    dest = backend(aggregate_model)
    src = backend(src_model)

    #MOI.empty!(dest)

    idxmap = MOI.Utilities.IndexMap()

    #Copy variables
    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    vars = MOI.add_variables(dest, length(vis_src))
    for (vi, var) in zip(vis_src, vars)
        idxmap.varmap[vi] = var
    end

    # Copy variable attributes
    MOI.Utilities.pass_attributes(dest, src, copy_names, idxmap, vis_src)

    # Copy model attributes
    MOI.Utilities.pass_attributes(dest, src, copy_names, idxmap)

    # Copy constraints
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        # do the rest in copyconstraints! which is type stable
        MOI.Utilities.copyconstraints!(dest, src, copy_names, idxmap, F, S)
    end

    return idxmap
end

    #jump_node.variablemap = node_map

#jump_node.indexmap = index_map

# for (key, data) in node_model.ext
#     new_model.ext[key] = copy_extension_data(data, new_model, model)
# end
#Setup new variable

#Upper and Lower Bounds
# JuMP.setlowerbound(x,JuMP.lower_bound(var))
# JuMP.setupperbound(x,JuMP.upper_bound(var))
#
                                 #rename the variable to the node model variable name plus the node or edge name
#
# #setcategory(x,node_model.colCat[i])                                  #set the variable to the same category
# if JuMP.is_binary(var)
#     JuMP.set_binary(new_x)
# end
#
# if JuMP.is_integar(var)
#     JuMP.set_integer(new_x)
# end
#
# if JuMP.is_fixed(var)
#     JuMP.set_fixed(var,) #Query fixed value
# end



#TODO
# function setsumgraphobjectives(graph)
#     has_nonlinear = false
#     #check if any nodes have nonlinear objectives
#     for node in getnodesandedges(graph)
#         node_model = getmodel(node)
#         traits = JuMP.ProblemTraints(node_model)
#         if traits.nlp == true
#             has_nonlinear = true
#             break
#         end
#     end
#     #if it's all linear or quadtratic
#     if has_nonlinear  == false
#         obj = 0
#         for node in getnodesandedges(graph)
#             node_model = getmodel(node)
#             obj += node_model.obj
#         end
#         graph.obj = obj  #set the graph objective to the sum of each node
#     elseif has_nonlinear == true
#         obj = :()
#         for node in getnodesandedges(graph)
#             node_model = getmodel(node)
#             d = JuMP.NLPEvaluator(node_model)
#             MathProgBase.initialize(d,[:ExprGraph])
#             node_obj = MathProgBase.obj_expr(d)
#             obj = Expr(:call,:+,copy(obj),node_obj)
#             # ex1 = MathProgBase.obj_expr(d)
#             # ex2 = MathProgBase.obj_expr(d2)
#             # newexpr = Expr(:call, :+, copy(ex1), copy(ex2))
#             # JuMP.setNLobjective(m, P.objSense, newexpr)
#         end
#     end
# end

# for (key,var) in getnodevariablemap(node)
#     var2 = node2[key]
#     if isa(var,JuMP.JuMPArray) || isa(var,Array)# || isa(var,JuMP.Variable)
#         vals = JuMP.getvalue(var)  #get value of the
#         Plasmo.setarrayvalue(var2,vals)  #the dimensions have to line up for arrays
#     elseif isa(var,JuMP.JuMPDict) || isa(var,Dict)
#         Plasmo.setvalue(var,var2)
#     elseif isa(var,JuMP.Variable)
#         JuMP.setvalue(var2,JuMP.getvalue(var))
#     else
#         error("encountered a variable type not recognized")
#     end
# end

# #Containers
# #TODO Create an objdict with the reproduced JuMP container types for each node
# for key in keys(node_model.objDict)  #this contains both variable and constraint references
#     if isa(node_model.objDict[key],Union{JuMP.JuMPArray{AbstractJuMPScalar},Array{AbstractJuMPScalar}})     #if the JuMP variable is an array or a JuMPArray
#         vars = node_model.objDict[key]
#         isa(vars,JuMP.JuMPArray) ? vars = vars.innerArray : nothing
#         dims = JuMP.size(vars)
#         node_map[key] = Array{AbstractJuMPScalar}(dims)
#         for j = 1:length(vars)
#             var = vars[j]
#             node_map[key][j] = var_map[linearindex(var)]
#         end
#     #reproduce the same mapping in a dictionary
#     elseif isa(node_model.objDict[key],JuMP.JuMPDict)
#         tdict = node_model.objDict[key].tupledict  #get the tupledict
#         d_tmp = Dict()
#         for dkey in keys(tdict)
#             d_tmp[dkey] = var_map[linearindex(tdict[dkey])]
#         end
#         node_map[key] = d_tmp
#
#     elseif isa(node_model.objDict[key],JuMP.AbstractJuMPScalar) #else it's a single variable
#         node_map[key] = var_map[linearindex(node_model.objDict[key])]
#     # else #objDict also has contraints!
#     #     error("Did not recognize the type of a JuMP variable $(node_model.objDict[key])")
#     end
# end

#num_vars = MathProgBase.numvar(node_model)
#num_vars = JuMP.num_variables(node_model)

#var_map = Dict{MOI.Index,JuMP.VariableRef}()    #this dict will map linear index of the node model variables to the new model JuMP variables {node var index => flat model JuMP.Variable}
#node_map = Dict()             #nodemap. {varkey => [var1,var2,...]}
#index_map = Dict{MOI.Index,JuMP.VariableRef}()    #{var index in node => var index in flat model}

#add the node model variables to the new model
#Look at JuMP copy method for ideas here
#for i = 1:num_vars

#This might be unnecssary since I'll grab all of the variable data going through the MOI constraints

# for i = 1:length(node_model.linconstr)
#     con = node_model.linconstr[i]
#     #t = collect(linearterms(con.terms))  #This is broken in julia 0.5
#     t = []
#     for terms in linearterms(con.terms)
#         push!(t,terms)
#     end
#     reference = @constraint(m, con.lb <= sum(t[i][1]*var_map[linearindex(t[i][2])] for i = 1:length(t)) + con.terms.constant <= con.ub)
#     # push!(getattribute(nodeoredge,:NodeData).constraintlist,reference)
#     push!(jump_node.constraintlist,reference)
# end

# #COPY QUADRATIC CONSTRAINTS
# for i = 1:length(node_model.quadconstr)
#     con = node_model.quadconstr[i]
#     #collect the linear terms
#     t = []
#     for terms in linearterms(con.terms.aff)
#         push!(t,terms)
#     end
#     qcoeffs =con.terms.qcoeffs
#     qvars1 = con.terms.qvars1
#     qvars2 = con.terms.qvars2
#     #Might be a better way to do this
#     if con.sense == :(==)
#         reference = @constraint(m,sum(qcoeffs[i]*var_map[linearindex(qvars1[i])]*var_map[linearindex(qvars2[i])] for i = 1:length(qcoeffs)) +
#         sum(t[i][1]*var_map[linearindex(t[i][2])] for i = 1:length(t)) + con.terms.aff.constant == 0)
#     elseif con.sense == :(<=)
#         reference = @constraint(m,sum(qcoeffs[i]*var_map[linearindex(qvars1[i])]*var_map[linearindex(qvars2[i])] for i = 1:length(qcoeffs)) +
#         sum(t[i][1]*var_map[linearindex(t[i][2])] for i = 1:length(t)) + con.terms.aff.constant <= 0)
#     elseif con.sense == :(>=)
#         reference = @constraint(m,sum(qcoeffs[i]*var_map[linearindex(qvars1[i])]*var_map[linearindex(qvars2[i])] for i = 1:length(qcoeffs)) +
#         sum(t[i][1]*var_map[linearindex(t[i][2])] for i = 1:length(t)) + con.terms.aff.constant >= 0)
#     end
#     # push!(getattribute(nodeoredge,:NodeData).constraintlist,reference)
#     push!(jump_node.constraintlist,reference)
# end

#define some setvalue functions for convenience when dealing with JuMP JuMPArray type
#dimension of jarr2 must be greater than jarr1
# function setarrayvalue(jarr1::JuMP.JuMPArray,jarr2::JuMP.JuMPArray)# = setvalue(jarr1.innerArray,jarr2.innerArray)
#     for i = 1:length(jarr2.innerArray)
#         JuMP.setvalue(jarr1[i],jarr2[i])
#     end
# end
#
# function setarrayvalue(jarr1::Array,jarr2::Array)# = setvalue(jarr1.innerArray,jarr2.innerArray)
#     for i = 1:length(jarr2)
#         JuMP.setvalue(jarr1[i],jarr2[i])
#     end
# end
#
# function setarrayvalue(jarr1::JuMP.JuMPArray,jarr2::Array) # = setvalue(jarr1.innerArray,jarr2)
#     for i = 1:length(jarr2)
#         JuMP.setvalue(jarr1.innerArray[i],jarr2[i])
#     end
# end
#
# function setarrayvalue(jarr1::Array,jarr2::JuMP.JuMPArray) # = setvalue(jarr1,jarr2.innerArray)
#     for i = 1:length(jarr2)
#         JuMP.setvalue(jarr1[i],jarr2.innerArray[i])
#     end
# end
#
# #define some getvalue and setvalue functions for dealing with JuMPDict objects.
# function setvalue(dict::Dict,jdict::JuMP.JuMPDict)
#     for key in keys(dict)
#         jdict.tupledict[key] = dict[key]
#     end
# end
#
# function setvalue(jdict::JuMP.JuMPDict,dict::Dict)
#     for key in keys(jdict.tupledict)
#         dict[key] = jdict.tupledict[key]
#     end
# end

# for var in vars
#     model_node = getnode(var)
#     var_index = JuMP.linearindex(var)                   #index in modelgraph node
#     node_index = getindex(model_graph,model_node)       #node index in modelgraph
#     jump_node = getnode(jump_graph,node_index)          #the corresponding jumpgraph jumpnode
#     flat_indexmap = jump_node.indexmap
#     indexmap[var] = flat_indexmap[var_index]            #modelnode variable => jumpnode variable index
# end
# t = []
# for terms in linearterms(linkconstraint.terms)
#     push!(t,terms)
# end

#Create a linear constraint from the LinkConstraints
#new_constraint = JuMP.ScalarConstraint()
#JuMP.addconstraint(jump_model,linkconstraint.func,linkconstraint.set)
#con_reference = @constraint(jump_model, linkconstraint.lb <= sum(t[i][1]*JuMP.Variable(jump_model,indexmap[(t[i][2])]) for i = 1:length(t)) + linkconstraint.terms.constant <= linkconstraint.ub)
#push!(jump_graph.linkconstraints,con_reference)
