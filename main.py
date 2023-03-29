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

TIME = 0

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

    for h in range(n):
        # Find a pendant pair
        [t, u] = pendentpair(Fnew, inew)
        if t is None and u is None:
            return True
        # This gives a candidate solution
        A.append(S[u - 1].copy())
        isu = Fnew({u})
        if isu == -464647101:
            return True
        s.append(isu)
        if s[-1] < minimum:
            minimum = s[-1]
            position_of_min = len(s) - 1
        S[t - 1] += S[u - 1]
        del inew[u]
        for x in range(len(S[u - 1])):
            S[u - 1][x] *= -1
    return False


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

            ax = F(Wi)
            bx = F({V[j]})
            if ax == -464647101 or bx == -464647101:
                return [None, None]
            keys[counter] = ax - bx
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


def f_wrapper(G, result, keep, formula, N):
    nodes_set = set(G.vs.indices)

    iters = 1
    reso = 0
    objf = "CPM"

    bv = abs(ig.community._community_leiden(G, objective_function=objf, n_iterations=iters,
                                            resolution=reso).modularity * len(nodes_set))
    def chkp(r):
        x = [keep[int(r.vs[i]["_nx_name"])] for i in r.vs.indices]
        if len(set(x)) >= N:
            res = set()
            for v, _ in Counter(x).most_common():
                if v * -1 not in res:
                    res.add(v)
            if len(res) == N and checksat(formula, set(res)):
                return True

            res = set()
            for v, _ in reversed(Counter(x).most_common()):
                if v * -1 not in res:
                    res.add(v)
            if len(res) == N and checksat(formula, set(res)):
                return True

        return False

    def f(s):

        s = tuple(sorted(s))
        s1 = tuple(nodes_set - set(s))

        if s in result or s1 in result:
            return bv - result[s] + result[s1]
        else:
            global TIME
            TIME += 1

        g = G.subgraph(s1)
        h = G.subgraph(s)

        gvc, gec, hvc, hec = g.vcount() ,g.ecount(), h.vcount(), h.ecount()

        if gvc < 1 or gec < 1 or hvc < 1 or hec < 1:
            return 0

        if (gvc >= N and chkp(g)) or (hvc >= N and chkp(h)):
            return -464647101

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


def solve(g, keep, f, N):
    result = {}

    G = ig.Graph.from_networkx(g)
    q = queyranne(f_wrapper(G, result, keep, f, N), list(range(0, g.number_of_nodes())))
    return q == True


def here(g, keep, f, N):
    return solve(g, keep, f, N)


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
    N = 9
    M = 3

    yes = 0
    TM = set()
    for ui in range(0, 1):
        #f = list(gen_unsat(N, int(N * M)))
        f = gen_sat(N, int(N * M), random_cnf(N))
        g, k = sat_to_clique(f)
        p = 0
        yes += here(g, k, f, N)
        TM.add(TIME)
        TIME = 0

    print(yes)
    print(min(TM), max(TM))

    # 2 3098
