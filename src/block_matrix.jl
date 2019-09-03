####################################
# Analysis Functions
####################################
#Create an incidence matrix representing the underlying hypergraph
function getincidencematrix(graph::ModelGraph;include_master_node = false)
    if include_master_node
        #Add master node entry and create simple edges for link variables
        error("master node in incidence matrix not yet supported")
    else
        return sparse(graph.hypergraph)
    end
end

#Create a sparse matrix representing the ModelGraph structure.
function getblockdata(graph::ModelGraph)

    #TODO NestedBlockMatrix
    if has_subgraphs(graph)
        return getnestedblockdata(graph)
    end

    n_link_constraints = num_linkconstraints(graph)

    master = getmastermodel(graph)
    n_master_constraints = ModelGraphs.num_all_constraints(master)

    n_rows_top_block = n_link_constraints + n_master_constraints #master constraints + link constraints
    n_rows_bottom_block = sum(ModelGraphs.num_node_constraints.(getnodes(graph)))

    n_rows = n_rows_top_block + n_rows_bottom_block #total number of block matrix rows

    n_master_variables = JuMP.num_variables(master)
    n_cols = JuMP.num_variables(graph)

    #Initialze sparse matrix
    A = spzeros(n_rows,n_cols)
    top_offset = n_rows_top_block

    Amaster = ModelGraphs.jump_jacobian_structure(master)              #master problem jacobian structure
    Alink = ModelGraphs.link_jacobian_structure(graph)                 #link matrix jacobian structure

    for pair in Amaster
        A[pair[1],pair[2]] = 1
    end

    for pair in Alink
        A[pair[1] + n_master_constraints,pair[2] + n_master_variables] = 1
    end

    node_row_offsets = Vector{Int64}()
    node_col_offsets = Vector{Int64}()
    # node_col_ranges = Dict()
    # node_row_ranges = Dict()
    node_views = []

    node_row_offset = top_offset + 1
    node_col_offset = n_master_variables + 1

    for node in getnodes(graph)
        n_vars = JuMP.num_variables(getmodel(node))
        n_cons = ModelGraphs.num_node_constraints(node)

        Anode = ModelGraphs.jump_jacobian_structure(getmodel(node))

        #get node indices
        start_row = node_row_offset
        end_row = start_row + n_cons - 1
        start_col = node_col_offset
        end_col = start_col + n_vars - 1

        Aview = view(A,start_row:end_row,start_col:end_col)
        push!(node_views,Aview)

        for pair in Anode
            Aview[pair[1],pair[2]] = 1
        end

        for (local_var,master_var) in node.linkvariablemap
            master_column = master_var.vref.index.value

            local_column = local_var.index.value
            local_constraint_indices = Aview[:,local_column]

            for local_index in local_constraint_indices.nzind
                A[start_row + local_index - 1,master_column] = 1
            end

        end

        # node_row_ranges[node] = start_row:end_row
        # node_col_ranges[node] = start_col:end_col


        node_row_offset += n_cons
        node_col_offset += n_vars
    end

    return A,node_views#node_row_ranges,node_col_ranges
end

getblockmatrix(graph::ModelGraph) = getblockdata(graph)[1]

function getdetailedblockmatrix(graph::ModelGraph)
end

function getdetailednestedblockmatrix(graph::ModelGraph)
end

#TODO: This is a little more challenging
function getnestedblockmatrix(graph::ModelGraph)
end
