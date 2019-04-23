
struct ProjectionMap
    node_map::Dict{Int64,Int64}
    edge_map::Dict{LightGraphs.AbstractEdge,LightGraphs.AbstractEdge}
end
ProjectionMap() = ProjectionMap(Dict{Int64,Int64}(),Dict{LightGraphs.AbstractEdge,LightGraphs.AbstractEdge}())

function Base.getindex(projection_map::ProjectionMap,node_index::Int64)
    return projection_map.node_map[node_index]
end
function Base.getindex(projection_map::ProjectionMap, edge_index::LightGraphs.AbstractEdge)
    return projection_map.edge_map[edge_index]
end
Base.broadcastable(projection_map::ProjectionMap) = Ref(projection_map)

function Base.setindex!(projection_map::ProjectionMap,node_index_1::Int64,node_index_2::Int64)
    projection_map.node_map[node_index_2] = node_index_1
end
function Base.setindex!(projection_map::ProjectionMap, edge_index_1::LightGraphs.AbstractEdge,edge_index_2::LightGraphs.AbstractEdge)
    projection_map.edge_map[edge_index_2] = edge_index_1
end

function Base.merge!(proj_map1::ProjectionMap,proj_map2::ProjectionMap)
    for (k,v) in proj_map2.node_map
        proj_map1.node_map_map[k] = v
    end
    for (k,v) in proj_map2.edge_map
        proj_map1.edge_map[k] = v
    end
end

#A Unipartite graph (Standard Graph) where nodes correspond to model nodes and edges correspond to links between model nodes.
#This structure can be used to convert a ModelGraph structure to a Shared Constraint structure
mutable struct NodeUnipartiteGraph <: AbstractModelGraph
    structuregraph::StructureGraphs.StructureGraph
    v_weights::Dict{Int64,Int64}                     #vertex weights
    e_weights::Dict{LightGraphs.AbstractEdge,Int64}  #edge weights
end

NodeUnipartiteGraph() = NodeUnipartiteGraph(StructureGraph(LightGraphs.Graph),Dict{StructureNode,Int64}(),Dict{StructureNode,Int64}())

StructureGraphs.create_node(graph::NodeUnipartiteGraph) = StructureNode()
StructureGraphs.create_edge(graph::NodeUnipartiteGraph) = StructureEdge()

function string(graph::NodeUnipartiteGraph)
    "Unipartite Graph\ngraph_id: "*string(getlabel(graph))*"\nnodes:"*string((length(getnodes(graph))))
end

#Convert Hypergraph ==> NodeUnipartite Graph
function NodeUnipartiteGraph(graph::ModelGraph)
    ugraph = NodeUnipartiteGraph()
    projection_map = ProjectionMap()

    #Add the model nodes to the Unipartite graph
    for node in getnodes(graph)
        idx = getindex(graph,node)
        new_node = create_node(ugraph)

        n_vars = JuMP.num_variables(getmodel(node))
        #n_vars = length(getmodel(node).colVal)

        add_node!(ugraph,new_node,index = idx)   #index should be the same as the original graph
        new_index = getindex(ugraph,new_node)
        ugraph.v_weights[new_index] = n_vars  #node weights are number of variables

        projection_map[new_index] = idx  #Map new node index to original node index
    end

    #Add the edges between nodes
    for edge in getedges(graph)
        hyperedge = getindex(graph,edge)
        vertices = hyperedge.vertices
        for i = 1:length(vertices)
            node_from = getnode(ugraph,vertices[i])
            other_vertices = vertices[i+1:end]
            for j = 1:length(other_vertices)
                node_to = getnode(ugraph,other_vertices[j])
                new_edge = add_edge!(ugraph,node_from,node_to)
                new_index = getindex(ugraph,new_edge)
                if !haskey(ugraph.e_weights,new_index)
                    ugraph.e_weights[new_index] = 1
                else
                    ugraph.e_weights[new_index] += length(edge.linkconstraints)  #edge weights are number of link constraints
                end

                projection_map[new_index] = hyperedge   #Map new simple edge to original hyperedge

            end
        end
    end

    #Add edges from subgraphs to the unipartite graph
    for subgraph in getsubgraphlist(graph)
        for edge in getedges(subgraph)
            hyperedge = getindex(subgraph,edge)
            vertices = hyperedge.vertices       #these are indices in the subgraph

            #NOTE Need indices in graph, not the subgraph
            subgraph_nodes = [getnode(subgraph,i) for i in vertices]
            graph_vertices = [getindex(graph,node) for node in subgraph_nodes]  #these are the vertices in the original graph (and hence, the ugraph)
            for i = 1:length(graph_vertices)
                node_from = getnode(ugraph,graph_vertices[i]) #Actually, we need the right node index from the graph, note the ugraph
                other_vertices = graph_vertices[i+1:end]
                for j = 1:length(other_vertices)
                    node_to = getnode(ugraph,other_vertices[j])
                    new_edge = add_edge!(ugraph,node_from,node_to)
                    new_index = getindex(ugraph,new_edge)
                    if !haskey(ugraph.e_weights,new_index)
                        ugraph.e_weights[new_index] = 1
                    else
                        ugraph.e_weights[new_index] += length(edge.linkconstraints)  #edge weights are number of link constraints
                    end

                    projection_map[new_index] = hyperedge

                end
            end
        end
    end

    return ugraph, projection_map
end



#A Unipartite graph (Standard Graph) where nodes correspond to link constraints and edges correspond to shared model nodes.
#This structure can be used to convert a ModelGraph structure to a Shared Model (Variable) structure
mutable struct LinkUnipartiteGraph <: AbstractModelGraph
    structuregraph::StructureGraphs.StructureGraph
    v_weights::Dict{Int64,Int64}                     #vertex weights
    e_weights::Dict{LightGraphs.AbstractEdge,Int64}  #edge weights
end

LinkUnipartiteGraph() = LinkUnipartiteGraph(StructureGraph(LightGraphs.Graph),Dict{StructureNode,Int64}(),Dict{StructureNode,Int64}())

StructureGraphs.create_node(graph::LinkUnipartiteGraph) = StructureNode()
StructureGraphs.create_edge(graph::LinkUnipartiteGraph) = StructureEdge()

function string(graph::LinkUnipartiteGraph)
    "Unipartite Graph\ngraph_id: "*string(getlabel(graph))*"\nnodes:"*string((length(getnodes(graph))))
end

function LinkUnipartiteGraph(graph::ModelGraph)
end

#A Bipartite graph (Standard Graph) where one set of nodes corresponds to
mutable struct ModelBipartiteGraph <: AbstractModelGraph
    structuregraph::StructureGraphs.StructureGraph
    part1::Vector{Int64}   #partition 1 of the bipartite graph
    part2::Vector{Int64}   #partition 2 of the bipartite graph
end

ModelBipartiteGraph() = BipartiteModelGraph(StructureGraph(LightGraphs.Graph()),Int64[],Int64[])

StructureGraphs.create_node(graph::ModelBipartiteGraph) = StructureNode()
StructureGraphs.create_edge(graph::ModelBipartiteGraph) = StructureEdge()

function string(graph::ModelBipartiteGraph)
    "Bipartite Graph\ngraph_id: "*string(getlabel(graph))*"\nnodes:"*string((length(getnodes(graph))))
end


#Convert Hypergraph ==> Bipartite Graph
function ModelBipartiteGraph(graph::ModelGraph)
    bgraph = BipartiteGraph()

    #model nodes
    for node in getnodes(graph)
        idx = getindex(graph,node)
        new_node = create_node(bgraph)
        add_node!(bgraph,new_node,index = idx)  #keep the same indices
        push!(bgraph.part1,idx)
    end

    #hyper edges
    for edge in getedges(graph)
        hyperedge = getindex(graph,edge)
        vertices = hyperedge.vertices
        constraint_node = add_node!(bgraph)
        push!(bgraph.part2,getindex(bgraph,constraint_node))
        #connect this "node" to the other nodes
        for vertex in vertices
            model_node = getnode(bgraph,vertex)
            add_edge!(bgraph,constraint_node,model_node)
        end
    end

    #TODO
    #subgraphs
    return bgraph
end


# #Functions that do graph transformations to facilitate different decomposition algorithms
# function convert_to_shared_variable(graph::ModelGraph)
# end
#
# function convert_to_shared_constraint(graph::ModelGraph)
# end

# #Convert JuMP Model ==> Bipartite Graph
# function getbipartitegraph(model::JuMP.Model)
#
# end
#
# #Convert JuMP Model ==> Unipartite Graph
# function getunipartitegraph(graph::JuMP.Model)
#
# end
