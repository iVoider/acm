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
        for n in v.neighbors():
            nxt[n.index] = 1.0
        transition.append(nxt)

    Q = np.array(transition)
    Q = Q/Q.sum(axis=1, keepdims=1)

    evals, evecs = np.linalg.eig(Q.T)
    evec1 = evecs[:, np.isclose(evals, 1)]
    # Since np.isclose will return an array, we've indexed with an array
    # so we still have our 2nd axis.  Get rid of it, since it's only size 1.
    evec1 = evec1[:, 0]
    stationary = evec1 / evec1.sum()

    # eigs finds complex eigenvalues and eigenvectors, so you'll want the real part.
    stationary = stationary.real
    return stationary


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
    N = 4
    M = 3

    yes = list()

    for ui in range(0, 100):
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
            print('shit')
            break
