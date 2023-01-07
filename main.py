import difflib
import itertools
import warnings

from sat import gen_unsat, gen_sat, random_cnf

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
    for h in range(n - 1):
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
    for i in range(n):
        vold = vnew
        Wi += [vold]
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


def f_wrapper(G, result, seed):
    nodes_set = set(G.vs.indices)

    def f(s):
        s = tuple(sorted(s))
        if s in result:
            return result[s]

        g = G.subgraph(nodes_set - set(s))

        if g.vcount() < 1 or g.ecount() < 1:
            return 0
        # for other types (like directed), use:
        val = ig.community._community_leiden(g).modularity
        result[s] = val
        return val

    return f


def solve(G, linear):
    def optimization(g, linear, seed):
        result = {}
        G = ig.Graph.from_networkx(g)
        q = queyranne(f_wrapper(G, result, seed), list(range(0, g.number_of_nodes())))
        ret = {}
        for key, value in q:
            ret[tuple(sorted(
                [(int(G.vs[i]["_nx_name"]) if not linear else tuple(G.vs[i]["_nx_name"])) for i in key]))] = value
        return ret

    return optimization(G, linear, None)


def sat_to_clique(formula):
    g = nx.Graph()
    map = {}

    node = 0

    for clause in formula:
        ins = (node, node + 1, node + 2)
        map[clause] = ins
        g.add_nodes_from(ins)
        node += 3

    while map:
        clause, nodes = map.popitem()
        for literal, node in zip(clause, nodes):
            for k, v in map.items():
                for i, j in zip(k, v):
                    if i * -1 != literal:
                        g.add_edge(node, j)
    return g

def testcase_linear(g):
    largest = len(g.nodes()) // 3

    a = solve(g, False)
    rm = 0
    doex = set(g.nodes()) - set([i[0] for i in a.keys()])
    a[(0,)] = a[(1,)]
    sortout = dict([((x,y), (a[(x,)] + a[(y,)])) for x,y in g.edges])

    c, _ = zip(*Counter.most_common(sortout))

    delu = 0

    for x,y in reversed(c):
        if g.has_edge(x, y):
         g.remove_edge(x, y)
        elif g.has_edge(y, x):
         g.remove_edge(y, x)

        check = {x, y}
        while check:
            n = check.pop()
            if n in g.nodes:
             if g.degree(n) < largest - 1:
                for nei in g.neighbors(n):
                    check.add(nei)
                delu += 1
                g.remove_node(n)

        rm += 1
        if ig.Graph.from_networkx(g).clique_number() < largest:
          break

    return rm, delu, doex

if __name__ == '__main__':
     while True:
      f = gen_sat(6, 6 * 4, random_cnf(6))
      #f = gen_unsat(4, 4 * 4)
      g = sat_to_clique(f)
      print(testcase_linear(g.copy()))

      # syort nodes by edges value
