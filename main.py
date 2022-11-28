import random
import time
import warnings

import networkx as nx

from generator import generate_problem
from graphcut import cut_edges
from maxcut import local_consistent_max_cut
from queyranne import queyranne

warnings.simplefilter(action='ignore', category=FutureWarning)

import itertools
import leidenalg
import igraph as ig


def f_wrapper(H, seed=None,
              partition_type=leidenalg.CPMVertexPartition):
    G = ig.Graph.from_networkx(H)
    nodes_set = set(G.vs.indices)

    def f(s):
        g = ig.Graph.subgraph(G, nodes_set - set(s))
        h = G.copy()
        h.delete_vertices(s)

        if g.vcount() < 1 or g.ecount() < 1:
            return 0
        return leidenalg.find_partition(g, partition_type, seed=seed).modularity + leidenalg.find_partition(h, partition_type, seed=seed).modularity

    return f


def solve(g, seed=None):
    q = queyranne(f_wrapper(g, seed), list(range(0, g.number_of_nodes())))
    C = q[:len(q)//2]
    C = [x[0] for x in C]
    k = nx.cut_size(g, C)
    return k


if __name__ == '__main__':
    G = nx.erdos_renyi_graph(n=20, p=0.2)
    H = G.copy()
    local_consistent_max_cut(G)
    print(cut_edges(G))
    print(solve(H, seed=random.getrandbits(16)))
