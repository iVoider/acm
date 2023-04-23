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

from collections import Counter


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


def prob(Gr):
    transition = []
    for v in Gr.vs:
        nxt = [0.0] * Gr.vcount()
        vl = 1.0 / len(v.neighbors())
        for n in v.neighbors():
            nxt[n.index] = vl
        transition.append(nxt)

    transition_matrix = np.array(transition)

    transition_matrix_transp = transition_matrix.T
    eigenvals, eigenvects = np.linalg.eig(transition_matrix_transp)
    close_to_1_idx = np.isclose(eigenvals, 1)
    target_eigenvect = eigenvects[:, close_to_1_idx]
    target_eigenvect = target_eigenvect[:, 0]
    stationary_distrib = target_eigenvect / sum(target_eigenvect)
    return stationary_distrib


def multipass(m_g):
    m_g = m_g.linegraph()
    r = list()
    p = prob(m_g)
    for v in m_g.vs.indices:
        map_ping = [0] * m_g.vcount()
        nb = m_g.vs[v].neighbors()
        for n in nb:
            map_ping[n.index] = 1
        m_g[v] = 1
        r.append(ig.VertexClustering(m_g, map_ping).modularity * p[v].real)
    return sorted(r)


if __name__ == '__main__':
    N = 6
    M = 3

    yes = list()

    for ui in range(0, 1):
        f = gen_unsat(N, N * M)

        # f = gen_sat(N, N * M, random_cnf(N))

        g, k = sat_to_clique(f)
        G = ig.Graph.from_networkx(g)

        h = g.copy()
        while True:
            double_edge_swap(h, max_tries=1000, nswap=1)
            if not nx.is_isomorphic(g, h):
                break

        H = ig.Graph.from_networkx(h)
        mapping = H.vs.indices
        random.shuffle(mapping)
        H = H.permute_vertices(mapping)

        D = ig.Graph.from_networkx(g)
        mapping = D.vs.indices
        random.shuffle(mapping)
        D = D.permute_vertices(mapping)

        jk = multipass(G)
        if jk == multipass(H) or jk != multipass(D):
            break
        print(jk)
