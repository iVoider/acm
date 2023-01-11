import difflib
import itertools
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
         #val = leidenalg.find_partition(g, partition_type=leidenalg.CPMVertexPartition, seed=random.getrandbits(16)).modularity
         val = ig.community._community_leiden(g, objective_function="CPM").modularity
        else:
         #val = leidenalg.find_partition(g, partition_type=leidenalg.ModularityVertexPartition, seed=random.getrandbits(16)).modularity
         val = ig.community._community_leiden(g, objective_function="modularity").modularity
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
            [(int(G.vs[i]["_nx_name"]) if not linear else tuple(G.vs[i]["_nx_name"])) for i in key]))] = value
    di = {}
    for e in g.edges:
        di[e] = (ret[(e[0],)] + ret[(e[1],)])
    return Counter.most_common(di)


def here(g, ksize, cpm = True):
    k = [(frozenset(x), y) for x, y in solve(g, False, cpm)]

    d = dict(k)

    # v = {}
    #
    # for key, value in sorted(d.items()):
    #     v.setdefault(value, []).append(key)
    #
    # print(v)

    z = 0
    for origin in k[:1000]:
        z += 1
        options = {origin[0]}
        cur = set()
        while options:
            next = options.pop()
            co = tuple(next)
            cur.add(co[0])
            cur.add(co[1])
            pos = frozenset(g.neighbors(co[0]))

            for c in cur:
                pos = pos & frozenset(g.neighbors(c))

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
        if len(cur) >= ksize and nx.density(g.subgraph(cur)) == 1:
            return True

    return False


if __name__ == '__main__':
 r = 0
 while True:
    f = gen_sat(15, 15 * 3, random_cnf(15))
    #f = gen_unsat(10, 10 * 5)
    g = sat_to_clique(f)
    t = 0
    while True:
      if here(g.copy(), 15 * 3, True) or here(g.copy(), 15 * 3, False):
          break
      else:
        t += 1
    r += 1
    print(r, t)
