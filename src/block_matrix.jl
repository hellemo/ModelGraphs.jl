####################################
# Analysis Functions
####################################
#Create an incidence matrix representing the underlying hypergraph
function getincidencematrix(graph::ModelGraph;include_master_node = false)
    if include_master_node
        error("master node in incidence matrix not yet supported")
    else
        return sparse(graph.hypergraph)
    end
end

#Create a sparse matrix representing the ModelGraph structure.
function getblockdata(graph::ModelGraph)

    if has_subgraphs(graph)
        return getnestedblockmatrix(graph)
    end

    n_link_constraints = num_linkconstraints(graph)

    master = getmastermodel(graph)
    n_master_constraints = ModelGraphs.num_all_constraints(master)

    n_rows_top_block = n_link_constraints + n_master_constraints #master constraints + link constraints
    n_rows_bottom_block = sum(ModelGraphs.num_node_constraints.(getnodes(graph)))

    n_rows = n_rows_top_block + n_rows_bottom_block

    n_master_variables = JuMP.num_variables(master)
    n_cols = JuMP.num_variables(graph)

    #Initialze sparse matrix
    A = spzeros(n_rows,n_cols)
    top_offset = n_rows_top_block

    #Fill the Master Block
    A[1:n_rows_top_block,1:n_master_variables] .= 1

    node_row_offsets = Vector{Int64}()
    node_col_offsets = Vector{Int64}()

    #Fill the diagonal and master column
    # push!(node_row_offsets,n_master_constraints)
    # push!(node_col_offsets,n_master_variables)
    node_col_ranges = Dict()
    node_row_ranges = Dict()

    node_row_offset = top_offset + 1
    node_col_offset = n_master_variables + 1

    for node in getnodes(graph)
        #index = getindex(graph,node)
        n_vars = JuMP.num_variables(getmodel(node))
        n_cons = ModelGraphs.num_node_constraints(node)

        #get node indices
        start_row = node_row_offset
        end_row = start_row + n_cons - 1
        start_col = node_col_offset
        end_col = start_col + n_vars - 1
        A[start_row:end_row,start_col:end_col] .= 1
        if !(isempty(node.linkvariablemap))
            A[start_row:end_row,1:n_master_variables] .= 1
        end

        node_col_ranges[node] = start_col:end_col
        node_row_ranges[node] = start_row:end_row

        node_row_offset += n_cons
        node_col_offset += n_vars
    end

    offset = n_master_constraints + 1
    for linkconstraint in getlinkconstraints(graph)
        nodes = getnodes(linkconstraint)
        for node in nodes
            col_range = node_col_ranges[node]
            A[offset,col_range] .= 1
        end
        offset += 1
    end

    return A,node_row_ranges,node_col_ranges
end

getblockmatrix(graph::ModelGraph) = getblockdata(graph)[1]

function getdetailedblockmatrix(graph::ModelGraph)
end

function getdetailednestedblockmatrix(graph::ModelGraph)
end

#TODO: This is a little more challenging
function getnestedblockmatrix(graph::ModelGraph)
end
