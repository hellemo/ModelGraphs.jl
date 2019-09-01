using Plots

function Plots.spy(graph::ModelGraph;markershape = :square,markersize = 20,colorbar = nothing, c = ColorGradient([:red,:blue]), kwargs...)
    block_matrix = getblockmatrix(graph)

    #Set values so we map diagonal entries to a different color
    # #TODO: Handle subgraphs (sub-blocks)
    # for i = 1:getnumnodes(graph)
    #    block_matrix[i+1,i+1] = 3   #diagonal blocks will be blue
    # end

    p = Plots.spy(block_matrix, markershape = markershape, markersize = markersize, colorbar = colorbar, c = c, kwargs...)
    n = getnumnodes(graph) + 1
    xlims!(p,(0,n + 1))
    ylims!(p,(0,n + 1))

    return p
end
