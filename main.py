import difflib
import itertools
import warnings
import time

import leidenalg

import numpy as np
import scipy.stats
from scipy.stats import mannwhitneyu, ttest_ind, zscore

from sat import gen_unsat, gen_sat, random_cnf, sat, solution, asol, mincore

warnings.simplefilter(action='ignore', category=FutureWarning)

import igraph as ig
import networkx as nx

from collections import OrderedDict, Counter


def here(G):
    candidates = {}

    for v in G.vs.indices:
        H = G.copy()
        H.delete_vertices(v)
        vcc = ig.VertexClustering(H, H.vs.indices)
        candidates[v] = vcc.modularity

    consider = {}

    for e in G.vs[min(candidates, key=candidates.get)].incident():
        H = G.copy()
        H.delete_edges(e)
        vcc = ig.VertexClustering(H, H.vs.indices)
        a = vcc.modularity

        H = G.copy()
        H.delete_vertices([e.source, e.target])
        vcc = ig.VertexClustering(H, H.vs.indices)
        b = vcc.modularity

        H = G.copy()
        mapping = H.vs.indices
        mapping[e.source] = e.target
        H.contract_vertices(mapping)
        vcc = ig.VertexClustering(H, H.vs.indices)
        cs = vcc.modularity

        consider[e.index] = (a, b, cs)

    G.delete_edges(min(consider, key=consider.get))

    return G


# check remove edge between two smallest nodes values

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


if __name__ == '__main__':
    N = 6
    size = 18

    start = time.time()

    yes = 0
    for ui in range(0, 1):
        f = gen_sat(N, size, random_cnf(N))
        # f = gen_unsat(N,size)
        g, k = sat_to_clique(f)
        G = ig.Graph.from_networkx(g)
        i = 0
        print(G.ecount(), G.vcount())
        while True:
            G = here(G)
            G.vs.select(_degree=0).delete()
            G.vs.select(_degree=size - 2).delete()
            if G.clique_number() < size:
                print(G.ecount(), G.vcount())
                break
            i += 1
        print(i)
