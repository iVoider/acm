import difflib
import itertools
import random
import warnings
import time

import leidenalg

import numpy as np
import scipy.stats
from networkx import double_edge_swap
from scipy.stats import mannwhitneyu, ttest_ind, zscore

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


def f_wrapper(G, result, vals):
    nodes_set = set(G.vs.indices)

    def f(s):

        s = tuple(sorted(s))

        if s in result:
            return result[s]

        g = G.subgraph(nodes_set - set(s))

        if g.vcount() < 1 or g.ecount() < 1:
            return 0
        # for other types (like directed), use:
        partition = 0

        for ve in g.vs:
            partition += vals[int(ve["_nx_name"])]

        result[s] = partition
        return partition

    return f

def here(g, size, keep):
    G = ig.Graph.from_networkx(g)

    vals = {}
    vals_c = {}

    for v in G.vs:
        H = G.copy()
        H.delete_vertices(v)
        vals[int(v["_nx_name"])] = ig.Graph.community_leiden(H, objective_function="CPM").modularity
        vals_c[int(v["_nx_name"])] = ig.Graph.community_leiden(H, objective_function="modularity").modularity



    result = {}

    q = queyranne(f_wrapper(G, result, vals), list(range(0, g.number_of_nodes())))

    result = {}

    q1 = queyranne(f_wrapper(G, result, vals_c), list(range(0, g.number_of_nodes())))

    return keep[int(G.vs[q[0][0][0]]["_nx_name"])], keep[int(G.vs[q1[0][0][0]]["_nx_name"])], keep[int(G.vs[q[-1][0][0]]["_nx_name"])],  keep[int(G.vs[q1[-1][0][0]]["_nx_name"])]

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
    size = 12

    yes = [0] * 4


    for ui in range(0, 100):
        f = gen_sat(N, size, random_cnf(N))
        #f = gen_unsat(N,size)
        g, k = sat_to_clique(f)
        a,b,c,d = here(g, size, k)

        yes[0] += asol(f, asum=[a]) != []
        yes[1] += asol(f, asum=[b]) != []
        yes[2] += asol(f, asum=[c]) != []
        yes[3] += asol(f, asum=[d]) != []



    print(yes)
