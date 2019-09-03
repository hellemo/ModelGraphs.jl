using Plots

#TODO: Handle subgraphs (sub-blocks)
function Plots.spy(graph::ModelGraph;markershape = :square,markersize = 1,colorbar = nothing, c = ColorGradient([:red,:blue]), kwargs...)

    A,node_row_range,node_col_range = getblockdata(graph)

    #Set values so we map diagonal entries to a different color
    #NOTE: Just get specific indices for each node
    # for node in getnodes(graph)
    #     row_range = node_row_range[node]
    #     col_range = node_col_range[node]
    #     A[row_range,col_range] .= 3   #diagonal blocks will be blue
    # end

    p = Plots.spy(A, markershape = markershape, markersize = markersize, colorbar = colorbar, c = c, kwargs...)

    n_cols = JuMP.num_variables(graph)

    n_link_constraints = num_linkconstraints(graph)

    master = getmastermodel(graph)
    n_master_constraints = ModelGraphs.num_all_constraints(master)

    n_rows_top_block = n_link_constraints + n_master_constraints #master constraints + link constraints
    n_rows_bottom_block = sum(ModelGraphs.num_node_constraints.(getnodes(graph)))

    n_rows = n_rows_top_block + n_rows_bottom_block

    xlims!(p,(0,n_cols + 1))
    ylims!(p,(0,n_rows + 1))

    return p
end
