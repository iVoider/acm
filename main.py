import random
import warnings
from enum import Enum

from generator import generate_problem

warnings.simplefilter(action='ignore', category=FutureWarning)

from networkx.algorithms.community import louvain_communities
from queyranne import QUEYRANNE
import networkx as nx
import itertools
import leidenalg
import igraph as ig

global_rand_seed = random.getrandbits(16)


# TODO: cache calculated partions for reuse
# Cross algo run, run in parallel use min max value fro one to another
def optimize(H, partition_type=leidenalg.CPMVertexPartition, lounvain=False):
    def f(V, s, params=None):
        if type(s) == int:
            s = [list(H.nodes)[s]]
        else:
            s = list(itertools.chain(*(i if isinstance(i, list) else (i,) for i in s)))
            s = [list(H.nodes)[n] for n in s]

        g = H.copy()
        g.remove_nodes_from(s)

        if g.number_of_nodes() < 1 or g.number_of_edges() < 1:
            return 0

        if lounvain:
            return 1.0 - nx.community.modularity(g, louvain_communities(g, seed=random.getrandbits(16)))
        else:
            ih = ig.Graph.from_networkx(g)
            part = leidenalg.find_partition(ih, partition_type, seed= global_rand_seed)
            # part = [[ih.vs[i]["_nx_name"] for i in j] for j in part]
            return 1.0 - ig.community._modularity(ih, part)

    return f


class QueyranneType(Enum):
    FIRST = 0
    SECOND = 1
    THIRD = 2


def solve(A, B, mapping, partition_type=leidenalg.CPMVertexPartition,
          lounvain=False, algo_type=QueyranneType.FIRST):
    def optimization(g):
        m = nx.adjacency_matrix(g).toarray()
        if algo_type == QueyranneType.FIRST:
            return QUEYRANNE(m, optimize(g, partition_type, lounvain))

    x, y = optimization(A), optimization(B)

    an = list(A.nodes)
    bn = list(B.nodes)

    res_x, val_x = zip(*x)
    res_y, val_y = zip(*y)

    inv_map = {v: k for k, v in mapping.items()}

    return set(zip([frozenset([an[i] for i in j]) for j in res_x], val_x)), set(zip([frozenset([inv_map[bn[i]] for i in j]) for j in res_y], val_y))


if __name__ == '__main__':
    N = 50
    D = 11

    while True:
     A, B, mapping = generate_problem(D, N, isomorphic=False)
     x_g = set()
     y_g = set()

     for n in range(0, 5):
      global_rand_seed = random.getrandbits(16)
      x,y = solve(A, B, mapping, leidenalg.ModularityVertexPartition, lounvain=False)
      x_g = x_g | x
      y_g = y_g | y

     q = set()
     l = set()
     for a,b in itertools.product(x_g,y_g):
        n1,v1 = a
        n2,v2 = b
        if v1 == v2:
         if n1 == n2:
          q.add(n1)
         elif len(n1) == len(n2):
          l.add((n1, n2))

     general = set()
     for e in q:
         for num in e:
             general.add(num)

     print(len(q), len(l), len(general))

