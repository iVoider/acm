import difflib
import itertools
import warnings
import leidenalg

import numpy as np
import scipy.stats
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


def f_wrapper(G, result, cpm=True):
    nodes_set = set(G.vs.indices)

    def f(s):

        s = tuple(sorted(s))

        if s in result:
            return result[s]

        g = G.subgraph(nodes_set - set(s))

        if g.vcount() < 1 or g.ecount() < 1:
            return 0
        # for other types (like directed), use:
        partition = ig.community._community_leiden(g, objective_function="CPM" if cpm else "modularity",
                                                   n_iterations=0).modularity

        result[s] = partition
        return partition

    return f


def solve(g, cpm=True):
    result = {}
    return queyranne(f_wrapper(g, result, cpm), list(range(0, g.vcount())))


def here(g, N, M):
    G = ig.Graph.from_networkx(g)
    candidates = {}

    for v in G.es.indices:
        H = G.copy()
        H.delete_edges([v])
        candidates[v] = ig.community._community_leiden(H, objective_function="CPM", n_iterations=0).modularity

    bestval = min(candidates.values())

    consider = {}
    for k, v in candidates.items():
        if v == bestval:
            H = G.copy()
            H.delete_edges([k])
            c = ig.community._community_leiden(H, objective_function="modularity", n_iterations=1).modularity * -1

            H = G.copy()
            H.delete_vertices([G.es[k].source])
            a = ig.community._community_leiden(H, objective_function="CPM", n_iterations=0).modularity
            a1 = ig.community._community_leiden(H, objective_function="modularity", n_iterations=1).modularity * -1

            H = G.copy()
            H.delete_vertices([G.es[k].target])
            b = ig.community._community_leiden(H, objective_function="CPM", n_iterations=0).modularity
            b1 = ig.community._community_leiden(H, objective_function="modularity", n_iterations=1).modularity * -1

            consider[k] = (a + b, c, abs(abs(a) - abs(b)), a1 + b1, abs(abs(a1) - abs(b1)))

    # G.delete_edges([min(consider, key=consider.get)])
    # res = (G.clique_number() >= N * M)
    return 1, len(consider.values()) - len(set(consider.values())), len(candidates.values()) - len(set(candidates.values()))


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
    M = 3

    yes = 0
    adup = set()
    zup = set()
    for ui in range(0, 10000):

        f = gen_sat(N, int(N * M), random_cnf(N))
        #f = gen_unsat(N,int(N * M))
        g, k = sat_to_clique(f)
        y, d, z = here(g, N, M)
        yes += y
        adup.add(z)
        zup.add(d)

    print(yes)
    print(sorted(adup))
    print(sorted(zup))

    #[462, 467, 470, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 510]
    #[0, 1, 2, 3, 4, 5, 6, 8, 9, 11, 12, 14, 16, 23, 36, 39, 61, 63]
