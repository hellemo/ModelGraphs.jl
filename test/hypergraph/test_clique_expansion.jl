# using Pkg
# Pkg.activate("..")

using NestedHyperGraphs
using LightGraphs

hyper = HyperGraph()

n1 = add_node!(hyper)
n2 = add_node!(hyper)
n3 = add_node!(hyper)
n4 = add_node!(hyper)

add_hyperedge!(hyper,1,2,3)
add_hyperedge!(hyper,3,4,1)
add_edge!(hyper,1,2)



hyper_sub1 = HyperGraph()
add_node!(hyper_sub1)
add_node!(hyper_sub1)
add_node!(hyper_sub1)
add_node!(hyper_sub1)
add_node!(hyper_sub1)
add_node!(hyper_sub1)

add_hyperedge!(hyper_sub1,1,2,3)
add_hyperedge!(hyper_sub1,1,2)
add_hyperedge!(hyper_sub1,4,1,3)


hyper_sub2 = HyperGraph()
add_node!(hyper_sub2)
add_node!(hyper_sub2)
add_node!(hyper_sub2)
add_node!(hyper_sub2)
add_node!(hyper_sub2)
add_node!(hyper_sub2)

add_hyperedge!(hyper_sub2,1,2,3)
add_hyperedge!(hyper_sub2,1,2,3,4,5,6)


add_subgraph!(hyper,hyper_sub1)
add_subgraph!(hyper,hyper_sub2)

sub_n1 = getnode(hyper_sub1,1)

sub_n2 = getnode(hyper_sub2,1)

add_hyperedge!(hyper,sub_n1,sub_n2)

clique_graph,projection_map = clique_expansion(hyper)
