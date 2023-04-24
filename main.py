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
import statistics

def f_wrapper(G, size, result):
    def f(V, s, params=None):

        if type(s) == int:
            s = [list(G.vs.indices)[s]]
        else:
            s = list(itertools.chain(*(i if isinstance(i, list) else (i,) for i in s)))
            s = [list(G.vs.indices)[n] for n in s]

        s = tuple(sorted(s))

        if s in result:
            return result[s]

        p = set(G.vs.indices) - set(s)

        transition = V.copy()

        for i in s:
            for j in p:
                transition[i, j] = 0.0

        mapping = [0] * G.vcount()
        for e in s:
            mapping[e] = 1

        result[s] = statistics.variance(list(prob(transition))) * ig.VertexClustering(G, mapping).modularity * -1
        return result[s]

    return f


def solve(G, size):
    result = {}

    transition = []

    for v in G.vs:
        nxt = [0.0] * G.vcount()
        for n in v.neighbors():
            nxt[n.index] = 1.0
        transition.append(nxt)

    Q = np.array(transition)
    io = QUEYRANNE(Q, f_wrapper(G, size, result))
    return io


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


def prob(transition):
    Q = np.array(transition)
    Q = Q / Q.sum(axis=1, keepdims=1)
    Q = np.nan_to_num(Q)
    evals, evecs = np.linalg.eig(Q.T)
    evec1 = evecs[:, np.isclose(evals, 1)]
    # Since np.isclose will return an array, we've indexed with an array
    # so we still have our 2nd axis.  Get rid of it, since it's only size 1.

    try:
     evec1 = evec1[:, 0]
    except:
        return [0.0] * transition.shape[0]
    stationary = evec1 / evec1.sum()

    # eigs finds complex eigenvalues and eigenvectors, so you'll want the real part.
    stationary = stationary.real
    return stationary


def multipass(m_g):
    r = list()
    m_g = m_g.linegraph()
    p = prob(m_g)
    for v in m_g.vs.indices:
        map_ping = [0] * m_g.vcount()
        nb = m_g.vs[v].neighbors()
        for n in nb:
            map_ping[n.index] = 1
        map_ping[v] = 2
        r.append(ig.VertexClustering(m_g, map_ping).modularity * p[v])
    return r


if __name__ == '__main__':
    N = 4
    M = 12

    yes = 0

    for ui in range(0, 100):
        f = gen_sat(N, M, random_cnf(N))
        g, k = sat_to_clique(f)
        G = ig.Graph.from_networkx(g)

        dc = dict.fromkeys(set(k.values()), 0)
        for s, v in solve(G, M):
            for e in s:
                dc[k[G.vs[e]["_nx_name"]]] += v * len(s)

        ry = max(dc, key = dc.get)
        yes += asol(f, asum=[ry]) != []
    print(yes)



