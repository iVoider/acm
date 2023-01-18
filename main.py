import difflib
import itertools
import operator
import warnings

import numpy as np
import sympy

from sat import gen_unsat, gen_sat, random_cnf, sat_to_clique, random_clause

warnings.simplefilter(action='ignore', category=FutureWarning)

import random
import leidenalg
import igraph as ig
import networkx as nx
from collections import OrderedDict, deque, defaultdict, Counter
import time


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


def f_wrapper(G, result, cpm):
    nodes_set = set(G.vs.indices)

    def f(s):
        s = tuple(sorted(s))
        if s in result:
            return result[s]

        g = G.subgraph(nodes_set - set(s))

        if g.vcount() < 1 or g.ecount() < 1:
            return 0
        if cpm:
         val = ig.community._community_leiden(g, objective_function="CPM", n_iterations=1).modularity
        else:
         val = ig.community._community_leiden(g, objective_function="modularity", n_iterations=1).modularity
        result[s] = val
        return val

    return f


def solve(g, linear, cpm = True):
    result = {}
    G = ig.Graph.from_networkx(g)
    q = queyranne(f_wrapper(G, result,  cpm), list(range(0, g.number_of_nodes())))
    ret = {}
    for key, value in q:
        ret[tuple(sorted(
            [(int(G.vs[i]["_nx_name"]) if not linear else tuple(G.vs[i]["_nx_name"])) for i in key]))[0]] = value

    di = {}

    for e in g.edges:
        di[e] = ret[e[0]] + ret[e[1]]

    return Counter.most_common(di), ret


def here(g, ksize, cpm = True):
    i,j = solve(g, False, cpm)
    k = [(frozenset(x), y) for x, y in i]

    d = dict(k)

    for origin in k:
        options = {origin[0]}
        cur = dict()

        pos = set(g.neighbors(tuple(origin[0])[0])) & set(g.neighbors(tuple(origin[0])[1]))

        while options:

            next = options.pop()
            co = tuple(next)
            if j[co[0]] > j[co[1]]:
             cur[co[0]] = None
             cur[co[1]] = None
            else:
             cur[co[1]] = None
             cur[co[0]] = None

            if len(pos) > 0:
                mx = -1
                mxv = -10e10
                for p in pos:
                    result = 0
                    for pc in cur:
                        result += d[frozenset({p, pc})]
                    if result > mxv:
                        mx = p
                        mxv = result
                for c in cur:
                    options.add(frozenset({mx, c}))

                pos = pos & set(g.neighbors(mx))

            if len(cur) == ksize:
               return True

    return False


if __name__ == '__main__':
 r = 0
 N = 10
 M = 4
 for p in range(0,10000):
    f = gen_sat(N, N * M, random_cnf(N))
    #f = gen_unsat(10, 10 * 5)
    g = sat_to_clique(f)
    t = 0
    while True:
      if here(g.copy(), N * M, True) or here(g.copy(), N * M, False):
          break
      else:
        t += 1
    r += 1
    print(r, t)


