import difflib
import itertools
import random
import warnings
import leidenalg

import numpy as np
from networkx import double_edge_swap, random_regular_graph

from q1 import QUEYRANNE
from sat import gen_unsat, gen_sat, random_cnf, sat, solution, asol, mincore

warnings.simplefilter(action='ignore', category=FutureWarning)

import igraph as ig
import networkx as nx

from collections import  Counter


def sat_to_clique(formula):
    g = nx.Graph()
    mapping = {}
    node = 0

    for clause in formula:
        for literal in clause:
            mapping[node] = literal
            node += 1

    nodes = list(zip(*[iter(range(0, 3 * len(formula)))] * 3))

    while nodes:
        clause = nodes.pop()
        for element in nodes:
            for x, y in itertools.product(clause, element):
                if mapping[x] != mapping[y] * -1:
                    g.add_edge(x, y)

    return g, mapping


def multipass(m_g):
    r = list()
    m_g = m_g.linegraph()
    for v in m_g.vs.indices:
        K = m_g.copy()
        sl = K.vcount()
        nxt = list()
        while True:
            map_ping = K.vs.indices
            nb = K.vs[v].neighbors()
            for n in nb:
                map_ping[n.index] = v
            K.contract_vertices(map_ping)
            map_ping = [0] * K.vcount()
            map_ping[v] = 1
            nxt.append(ig.VertexClustering(K, map_ping).modularity * len(nb))
            l = len(K.vs[v].neighbors())
            if l > sl or l < 2:
                break
            break

        r.append(tuple(nxt))
    return sorted(r)


if __name__ == '__main__':
    N = 4
    M = 7

    mx = set()

    yes = 0
    no = 0

    for ui in range(0, 1):
        #f = gen_unsat(N, N * M)
        f = gen_sat(N, N * M, random_cnf(N))
        g, k = sat_to_clique(f)
        h = g.copy()

        G = ig.Graph.from_networkx(g)

        # while True:
        #     double_edge_swap(h, max_tries=1000, nswap=1)
        #     if not nx.is_isomorphic(g, h):
        #         break
        #
        # H = ig.Graph.from_networkx(h)
        # mapping = H.vs.indices
        # random.shuffle(mapping)
        # H = H.permute_vertices(mapping)

        # D = ig.Graph.from_networkx(g)
        # mapping = D.vs.indices
        # random.shuffle(mapping)
        # D = D.permute_vertices(mapping)

        mg = multipass(G)
        # yes += mg == multipass(H)
        # no += mg == multipass(D)

        z = Counter(mg).most_common()[-1][1]
        mx.add(z)

    print(sorted(mx))
    print(yes, no)
