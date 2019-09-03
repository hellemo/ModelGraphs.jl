using Plots

#TODO: Handle subgraphs (sub-blocks)
function plotblockmatrix(graph::ModelGraph;markershape = :square,colorbar = nothing, c = ColorGradient([:red,:blue]), kwargs...)

    A,node_views = getblockdata(graph)


    for node_view in node_views
        node_view[node_view .== 1] .== 3
    end

    kwargs2 = Dict()
    kwargs = Dict(kwargs)
    for (symbol,value) in kwargs
        kwargs2[symbol] = value
    end
    kwargs2[:markershape] = markershape
    kwargs2[:colorbar] = colorbar
    kwargs2[:c] = c

    p = Plots.spy(A;kwargs2...)

    n_cols = JuMP.num_variables(graph)

    n_link_constraints = num_linkconstraints(graph)

    master = getmastermodel(graph)
    n_master_constraints = ModelGraphs.num_all_constraints(master)

    n_rows_top_block = n_link_constraints + n_master_constraints #master constraints + link constraints
    n_rows_bottom_block = sum(ModelGraphs.num_node_constraints.(getnodes(graph)))

    n_rows = n_rows_top_block + n_rows_bottom_block

    xlims!(p,(0,n_cols + 1))
    ylims!(p,(0,n_rows + 1))

    xlabel!(p,"# Variables")
    ylabel!(p,"# Constraints")

    return p
end
