import difflib
import itertools
import warnings

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
        val = leidenalg.find_partition(g, partition_type= leidenalg.CPMVertexPartition).modularity
        result[s] = val
        return val

    return f


def solve(G,H, linear):
    def optimization(g, linear, seed):
        result = {}
        G = ig.Graph.from_networkx(g)
        q = queyranne(f_wrapper(G, result, seed), list(range(0, g.number_of_nodes())))
        ret = {}
        for key, value in q:
            ret[tuple(sorted(
                [(int(G.vs[i]["_nx_name"]) if not linear else tuple(G.vs[i]["_nx_name"])) for i in key]))] = value
        return ret

    seed = random.getrandbits(16)

    x, y = optimization(G, False, seed), optimization(H, False, None)

    d = {n: [k for k in x.keys() if x[k] == n] for n in set(x.values())}

    print(d)

    d = {n: [k for k in y.keys() if y[k] == n] for n in set(y.values())}

    print(d)

    q = 0

    for at, bt in itertools.product(x.items(), y.items()):
        k1, v1 = at
        k2, v2 = bt
        if v1 == v2 and len(k1) == len(k2):
            q += 1

    return q

def generate_problem(degree, size, isomorphic):
    A = nx.random_regular_graph(degree, size)
    node_mapping = dict(zip(A.nodes(), sorted(A.nodes(), key=lambda k: random.random())))
    B = nx.relabel_nodes(A, node_mapping)

    k = list(B.nodes)
    random.shuffle(k)
    H = nx.Graph()
    H.add_nodes_from(k)
    H.add_edges_from(B.edges(data=True))
    B = H

    A.remove_nodes_from(list(nx.isolates(A)))
    B.remove_nodes_from(list(nx.isolates(B)))

    if not isomorphic:
        B = nx.double_edge_swap(B)
        # B = directed_edge_swap(B)
    if isomorphic != nx.is_isomorphic(A, B):
        return generate_problem(degree, size, isomorphic)

    return A, B, nx.isomorphism.vf2pp_isomorphism(A, B)

if __name__ == '__main__':
    G,H, m = generate_problem(5, 30, True)
    print(m)
    print(solve(G,H, False))
