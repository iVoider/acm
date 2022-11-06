# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


import warnings

import numpy as np
import scipy.linalg
from community import modularity
from networkx.algorithms.community import louvain_communities
from scipy.spatial.distance import directed_hausdorff
import sage as sp

from Queyanne import QUEYRANNE
from Sandpile import sandpile

warnings.simplefilter(action='ignore', category=FutureWarning)

import networkx as nx
import random
import itertools


def random_edge(nodes, ignore, used_edges, start):
    while True:
        choice = random.choice(nodes)
        if choice != ignore and (start, choice) not in used_edges:
            return start, choice


def generate_rand_digraphs(N, sparse=0.3, dropout=8, isomorphic=True):
    A = nx.MultiDiGraph()
    B = nx.MultiDiGraph()
    nodesA = list(range(0, N))
    nodesB = nodesA.copy()
    random.shuffle(nodesB)
    mapping = dict(zip(nodesA, nodesB))
    ins = dict(zip(nodesB, [set() for _ in range(N)]))
    outs = dict(zip(nodesB, [set() for _ in range(N)]))

    insert = []

    while len(A.nodes) < N:
        x = random.randint(0, N - 1)
        y = random.randint(0, N - 1)

        if random.random() < sparse:
            A.add_edge(x, y)
            insert.append((mapping[x], mapping[y]))
            ins[mapping[x]].add(mapping[y])
            outs[mapping[y]].add(mapping[x])

        if random.random() < sparse:
            A.add_edge(y, x)
            insert.append((mapping[y], mapping[x]))
            ins[mapping[y]].add(mapping[x])
            outs[mapping[x]].add(mapping[y])

    random.shuffle(insert)

    for edge in insert:
        B.add_edge(edge[0], edge[1])

    if not isomorphic:
        for x, y in itertools.combinations(nodesB, 2):
            if len(ins[x]) == len(ins[y]) and len(outs[x]) == len(outs[y]):
                if len(ins[x]) > 0 and len(ins[y]) > 0:
                    a = ins[x].pop()
                    b = ins[y].pop()
                    ins[x].add(b)
                    ins[y].add(a)
                    B.remove_edge(x, a)
                    B.add_edge(x, b)
                    B.remove_edge(y, b)
                    B.add_edge(y, a)
                    dropout -= 1
                elif len(outs[x]) > 0 and len(outs[y]) > 0:
                    a = outs[x].pop()
                    b = outs[y].pop()
                    outs[x].add(b)
                    outs[y].add(a)
                    B.remove_edge(a, x)
                    B.add_edge(b, x)
                    B.remove_edge(b, y)
                    B.add_edge(a, y)
                    dropout -= 1
            if dropout == 0:
                break

        if dropout != 0:
            return generate_rand_digraphs(N, sparse, dropout, isomorphic)

    return A, B, mapping


def m_cutfun(H):
    def f(V, s):
        if type(s) == int:
            s = [list(H.nodes)[s]]
        else:
            s = list(itertools.chain(*(i if isinstance(i, list) else (i,) for i in s)))
            s = [list(H.nodes)[n] for n in s]

        g = H.subgraph(s)
        if g.number_of_nodes() < 1 or g.number_of_edges() < 1:
            return 0
        return nx.algorithms.community.modularity(g, louvain_communities(g))

    return f


def solve(A, B, mapping):
    def unwind(g):
        m = nx.adjacency_matrix(g).toarray()
        cutfun = m_cutfun
        return QUEYRANNE(m, cutfun(g))

    x, y = unwind(A), unwind(B)

    an = list(A.nodes)
    bn = list(B.nodes)

    res_x, val_x = zip(*x)
    res_y, val_y = zip(*y)

    inv_map = {v: k for k, v in mapping.items()}

    return list(zip([[an[i] for i in j] for j in res_x], val_x)), list(
        zip([[inv_map[bn[i]] for i in j] for j in res_y], val_y))


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    A, B, mapping = generate_rand_digraphs(40, dropout=3, sparse=0.6, isomorphic=True)
    print(nx.is_isomorphic(A, B))
    x, y = solve(A, B, mapping)
    print(x)
    print(y)
