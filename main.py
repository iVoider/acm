import difflib
import itertools
import random
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


def here(g, k):
    G = ig.Graph.from_networkx(g)
    pem = list(range(0, len(G.vs.indices)))
    random.shuffle(pem)
    G = G.permute_vertices(pem)

    candidates = {}

    for v in G.vs.indices:
        H = G.copy()
        H.delete_vertices(v)
        candidates[v] = ig.community._community_leiden(H, objective_function="CPM",
                                                       n_iterations=1, resolution=0).modularity

    z = sorted(candidates.values())
    return abs(z[0] - z[-1]) * 10e15


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
    N = 5
    size = N * 3

    start = time.time()

    dd = list()

    for ui in range(0, 100):
        f = gen_sat(N, size, random_cnf(N))
        # f = gen_unsat(N,size)
        g, k = sat_to_clique(f)
        for z in range(0, 10):
         if here(g, k) != here(g, k):
            print("Shit T")
        dd.append((here(g, k), 1))

    for ui in range(0, 100):
        f = gen_unsat(N, size)
        g, k = sat_to_clique(f)
        for z in range(0, 10):
            if here(g, k) != here(g, k):
                print("Shit F")
        dd.append((here(g, k), 0))

    shit = {}
    shits = 0

    for x, y in dd:
        if x not in shit:
            shit[x] = {y}
        else:
            shit[x].add(y)
        if len(shit[x]) > 1:
            shits += 1

    dd = sorted(dd)
    dx, dc = zip(*dd)
    print(len(set(dx)), len(dx), set(dx), shits)

    vl = {}
    l = set(dx)
    vl[l.pop()] = [0, 0]
    vl[l.pop()] = [0, 0]
    for x, c in zip(dx, dc):
        vl[x][c] += 1
    print(vl)
    # import plotly.express as px
    # fig = px.scatter(x=dx, color=dc)
    # fig.show()
