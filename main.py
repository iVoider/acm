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

def here(g, size):
    G = ig.Graph.from_networkx(g)
    candidates = {}

    for v in G.vs.indices:
        H = G.copy()
        H.delete_vertices(v)
        candidates[v] = ig.community._community_leiden(H, objective_function="CPM", n_iterations=0).modularity / len(G.neighbors(v))

    G.delete_vertices(min(candidates, key=candidates.get))

    return G.clique_number() >= size

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
    N = 4
    size = 13

    start = time.time()

    yes = 0
    for ui in range(0, 100000):
        f = gen_sat(N, size, random_cnf(N))
        # f = gen_unsat(N,size)
        g, k = sat_to_clique(f)
        yes += here(g, size)

    print(yes)
    end = time.time()
    print(end - start)
