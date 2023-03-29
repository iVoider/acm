import difflib
import itertools
import warnings
import leidenalg

import numpy as np

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


def f_wrapper(G, result):
    nodes_set = set(G.vs.indices)

    iters = 1
    reso = 0
    objf = "CPM"

    bv = abs(ig.community._community_leiden(G, objective_function=objf, n_iterations=iters,
                                            resolution=reso).modularity * len(nodes_set))

    def f(s):

        s = tuple(sorted(s))
        s1 = tuple(nodes_set - set(s))

        if s in result or s1 in result:
            return bv - result[s] + result[s1]

        g = G.subgraph(s1)
        h = G.subgraph(s)

        if g.vcount() < 1 or g.ecount() < 1 or h.vcount() < 1 or h.ecount() < 1:
            return 0

        result[s1] = abs(ig.community._community_leiden(
            g,
            objective_function=objf,
            n_iterations=iters,
            resolution=reso).modularity * g.vcount())

        result[s] = abs(ig.community._community_leiden(
            h,
            objective_function=objf,
            n_iterations=iters,
            resolution=reso).modularity * h.vcount())

        partition = bv - result[s1] + result[s]
        result[s] = partition
        return partition

    return f


def solve(g, keep):
    result = {}

    G = ig.Graph.from_networkx(g)
    q = queyranne(f_wrapper(G, result), list(range(0, g.number_of_nodes())))
    ret = {}
    for key, value in result.items():
        to = [keep[int(G.vs[i]["_nx_name"])] for i in key]
        ret[tuple(to)] = value

    return ret


def here(g, keep):
    k = solve(g, keep)

    return Counter.most_common(k)


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


def checksat(f, ass):
    for c in f:
        if len(ass & set(c)) == 0:
            return False
    return True


if __name__ == '__main__':
    N = 4
    M = 3

    yes = 0
    mx = set()
    for ui in range(0, 100):
        # f = list(gen_unsat(N, int(N * M)))
        f = gen_sat(N, int(N * M), random_cnf(N))
        g, k = sat_to_clique(f)
        p = 0
        ans = here(g, k)
        mx.add(len(ans))
        for j, v in ans:
            p += 1
            if len(set(j)) >= N:
                res = set()
                for v, _ in reversed(Counter(j).most_common()):
                    if v * -1 not in res:
                        res.add(v)
                if len(res) == N and checksat(f, set(res)):
                    yes += 1
                    break

    print(yes)
    print(min(mx), max(mx))

    #         if Counter([int(g.vs[i]["_nx_name"]) for i in s1]).most_common():
