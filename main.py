import random
import warnings
from enum import Enum

warnings.simplefilter(action='ignore', category=FutureWarning)

from generator import generate_rand_digraphs
from networkx.algorithms.community import louvain_communities
from queyranne import QUEYRANNE
import networkx as nx
import itertools
import leidenalg
import igraph as ig


def optimize(H, partition_type=leidenalg.CPMVertexPartition, lounvain=False):
    def f(V, s, params=None):
        if type(s) == int:
            s = [list(H.nodes)[s]]
        else:
            s = list(itertools.chain(*(i if isinstance(i, list) else (i,) for i in s)))
            s = [list(H.nodes)[n] for n in s]

        g = H.subgraph(s)
        if g.number_of_nodes() < 1 or g.number_of_edges() < 1:
            return 0

        if lounvain:
            return 1.0 - nx.community.modularity(g, louvain_communities(g))
        else:
            ih = ig.Graph.from_networkx(g)
            part = leidenalg.find_partition(ih, partition_type)
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

    return list(zip([[an[i] for i in j] for j in res_x], val_x)), list(
        zip([[inv_map[bn[i]] for i in j] for j in res_y], val_y))


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    for n in range(0, 1):
        A, B, mapping = generate_rand_digraphs(20, dropout=random.randint(1, 3),
                                               sparse=random.uniform(0.2, 0.95),
                                               isomorphic=True)
        x, y = solve(A, B, mapping, leidenalg.ModularityVertexPartition)
        print(x)
        print(y)

        x, y = solve(A, B, mapping, leidenalg.ModularityVertexPartition)
        print(x)
        print(y)
