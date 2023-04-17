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

from collections import OrderedDict, Counter


def queyranne(F, V):
    def Fnew(a):
        r = []
        for x in a:
            r += S[x - 1]
        return F(r)

    n = len(V)
    S = [[x] for x in V]
    s = []
    A = []
    inew = OrderedDict()
    for x in range(1, n + 1):
        inew[x] = x
    minimum = float("inf")
    position_of_min = 0

    for h in range(n):
        # Find a pendant pair
        [t, u] = pendentpair(Fnew, inew)
        # This gives a candidate solution
        A.append(S[u - 1].copy())
        s.append(Fnew({u}))
        if s[-1] < minimum:
            minimum = s[-1]
            position_of_min = len(s) - 1
        S[t - 1] += S[u - 1]
        del inew[u]
        for x in range(len(S[u - 1])):
            S[u - 1][x] *= -1
    vals = dict(zip([tuple(a) for a in A], s))
    return Counter.most_common(vals)


# Implements the pendant pair finding subroutine of Queyranne's algorithm
# (Queyranne '95)
# F is the submodular function
# inds is an array of indices; (typically, 1:n)

def pendentpair(F, V):
    vstart = V.popitem(last=False)[0]
    vnew = vstart
    n = len(V)
    Wi = []
    used = [0] * n
    for i in range(n + 1):
        vold = vnew
        Wi += [vold]
        # Now update the keys
        keys = [1e99] * n
        minimum = float("inf")
        counter = -1
        for j in V:
            counter += 1
            if used[counter]:
                continue
            Wi += [V[j]]
            keys[counter] = F(Wi) - F({V[j]})
            del Wi[-1]
            if keys[counter] < minimum:
                minimum = keys[counter]
                argmin_key = j
                argmin_position = counter
            vnew = argmin_key
            used[argmin_position] = 1
    V[vstart] = vstart
    V.move_to_end(vstart, last=False)
    return [vold, vnew]


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

        vex = [0] * (size * 3)

        for i in s:
            vex[i] = 1

        result[s] = 1.0 - ig.VertexClustering(G, membership=vex).modularity
        return result[s]

    return f


def solve(G, size, g):
    result = {}
    return QUEYRANNE(nx.adjacency_matrix(g).toarray(), f_wrapper(G, size, result))


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
    M = 3

    yes = 0

    for ui in range(0, 100):
        #f = gen_unsat(N, N * M)
        f = gen_sat(N, N * M, random_cnf(N))
        g, k = sat_to_clique(f)

        node_mapping = dict(zip(g.nodes(), sorted(g.nodes(), key=lambda k: random.random())))
        h = nx.relabel_nodes(g, node_mapping)
        #
        # while True:
        #     double_edge_swap(h, max_tries=1000, nswap=1)
        #     if not nx.is_isomorphic(g, h):
        #         break

        G = ig.Graph.from_networkx(g)
        H = ig.Graph.from_networkx(h)

        yes += solve(G, N * M, g) == solve(H, N * M, h)


    print(yes)
