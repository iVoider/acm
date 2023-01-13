import difflib
import itertools
import operator
import warnings

from sat import gen_unsat, gen_sat, random_cnf, sat_to_clique

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
    for x in range(0, n + 1):
        inew[x] = x
    minimum = (float("inf"),float("inf"))
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


def pendentpair(F, V):
    vstart = V.popitem(last=False)[0]
    vnew = vstart
    n = len(V)
    Wi = []
    used = [0] * n
    for i in range(n + 1):
        vold = vnew
        Wi += [vold]
        keys = [1e99] * n
        minimum = (float("inf"), float("inf"))
        counter = -1
        for j in V:
            counter += 1
            if used[counter]:
                continue
            Wi += [V[j]]
            keys[counter] = tuple(map(operator.add, F(Wi), F({V[j]})))
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


def f_wrapper(G, result, cpm = True):
    nodes_set = set(G.vs.indices)

    def f(s):
        s = tuple(sorted(s))
        if s in result:
            return result[s]

        g = G.subgraph(nodes_set - set(s))

        if g.vcount() < 1 or g.ecount() < 1:
            return 0
        # for other types (like directed), use:
        if cpm:
         val = (ig.community._community_leiden(g, objective_function="CPM").modularity, 0.0)
        else:
         val = (ig.community._community_leiden(g, objective_function="modularity").modularity, 0.0)
        result[s] = val
        return val

    return f


def solve(g, linear, cpm = True):
    result = {}
    G = ig.Graph.from_networkx(g)
    q = queyranne(f_wrapper(G, result, cpm), list(range(0, g.number_of_nodes())))
    ret = {}
    for key, value in q:
        ret[tuple(sorted(
            [(int(G.vs[i]["_nx_name"]) if not linear else tuple(G.vs[i]["_nx_name"])) for i in key]))[0]] = value
    return ret


def here(g, ksize, cpm = True):
    x = solve(g, False, cpm)
    y = solve(g, False, not cpm)

    a = min(x, key=x.get)
    b = max(y, key=y.get)

    h = g.copy()
    h.remove_node(a)

    v = g.copy()
    v.remove_node(b)

    return ig.Graph.from_networkx(h).clique_number() == ksize or ig.Graph.from_networkx(v).clique_number() == ksize

if __name__ == '__main__':

    #f = gen_unsat(10, 10 * 5)

    y = 0
    for i in range(0, 10000):
     f = gen_sat(5, 5 * 4, random_cnf(4))
     g = sat_to_clique(f)
     if here(g, 5 * 4):
         y += 1
    print(y)



