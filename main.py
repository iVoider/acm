import random
import time
import warnings
from enum import Enum

from generator import generate_problem
from queyranne2 import queyranne

warnings.simplefilter(action='ignore', category=FutureWarning)

import itertools
import leidenalg
import igraph as ig


def f_wrapper(H, seed=None,
              partition_type=leidenalg.ModularityVertexPartition):
    nodes = list(H.nodes)

    def f(s):

        g = H.copy()
        g.remove_nodes_from(s)

        if g.number_of_nodes() < 1 or g.number_of_edges() < 1:
            return 0

        ih = ig.Graph.from_networkx(g)
        part = leidenalg.find_partition(ih, partition_type, seed=seed)
        return 1.0 - ig.community._modularity(ih, part)

    return f


def solve(A, B, seed=None):
    def optimization(g):
        return queyranne(f_wrapper(g, seed), list(range(0, g.number_of_nodes())))

    x, y = optimization(A), optimization(B)

    q = 0

    for at, bt in itertools.product(x, y):
        k1, v1 = at
        k2, v2 = bt
        if v1 == v2 and len(k1) == len(k2):
            q += 1

    return q


if __name__ == '__main__':
    N = 50
    D = 5

    A, B, mapping = generate_problem(D, N, isomorphic=False, directed=True)
    print(solve(A, B, seed=random.getrandbits(16)))
