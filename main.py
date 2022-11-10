import random
import time
import warnings

from generator import generate_problem
from queyranne import queyranne

warnings.simplefilter(action='ignore', category=FutureWarning)

import itertools
import leidenalg
import igraph as ig


def f_wrapper(H, seed=None,
              partition_type=leidenalg.ModularityVertexPartition):

    G = ig.Graph.from_networkx(H)
    nodes_set = set(G.vs.indices)
    def f(s):

        g = ig.Graph.subgraph(G, nodes_set - set(s))

        if g.vcount() < 1 or g.ecount() < 1:
            return 0

        return 1.0 - leidenalg.find_partition(g, partition_type, seed=seed).modularity

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
    N = 100
    D = 9

    A, B, mapping = generate_problem(D, N, isomorphic=True, directed=True)

    start = time.time()
    print(solve(A, B, seed=random.getrandbits(16)))
    print(time.time() - start)
