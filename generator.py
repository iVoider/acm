from random import random

import networkx as nx
from networkx import random_regular_graph, double_edge_swap

def generate_problem(degree, size, isomorphic):
    A = random_regular_graph(degree, size)
    node_mapping = dict(zip(A.nodes(), sorted(A.nodes(), key=lambda k: random())))
    B = nx.relabel_nodes(A, node_mapping)

    A.remove_nodes_from(list(nx.isolates(A)))
    B.remove_nodes_from(list(nx.isolates(B)))

    if not isomorphic:
        B = double_edge_swap(B)
    if isomorphic != nx.is_isomorphic(A, B):
        return generate_problem(degree, size, isomorphic)

    return A, B, node_mapping
