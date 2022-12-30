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

        g = ig.Graph.subgraph(G, nodes_set - set(s))

        if g.vcount() < 1 or g.ecount() < 1:
            return 0
        # for other types (like directed), use:
        val = 1.0 - leidenalg.find_partition(g, partition_type=leidenalg.ModularityVertexPartition,
                                             seed=seed).modularity
        result[s] = val
        return val

    return f


def are_isomorphic(G, H, linear):
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
    A = optimization(G, linear, seed)
    B = optimization(H, linear, seed)
    return A, B


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


def generate_problem_sub(degree, size, isomorphic):
    while True:
        try:
            A = nx.random_regular_graph(degree, size)
        except:
            A = None
        if A != None:
            break

    node_mapping = dict(zip(A.nodes(), sorted(A.nodes(), key=lambda k: random.random())))
    B = nx.relabel_nodes(A, node_mapping)

    rm = random.sample(list(B.nodes()), 10)

    B.remove_nodes_from(rm)

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
    matcher = nx.isomorphism.GraphMatcher(A, B)
    if isomorphic != matcher.subgraph_is_isomorphic():
        return generate_problem_sub(degree, size, isomorphic)
    return A, B, matcher.mapping


def topological_nodes(G, GN):
    VN = dict.fromkeys(G.nodes, (0, 0))
    for k, v in GN.items():
        for e in k:
            VN[e] = (VN[e][0] + v, VN[e][1] + 1)
    for k in VN:
        VN[k] = VN[k][0] / VN[k][1]
    a,b = zip(*Counter.most_common(VN))
    return a


def topological_edges(G, GE):
    VN = dict.fromkeys(G.edges, (0, 0))
    for k, v in GE.items():
        for e in k:
            VN[e] = (VN[e][0] + v, VN[e][1] + 1)
    for k in VN:
        VN[k] = VN[k][0] / VN[k][1]
    a, b = zip(*Counter.most_common(VN))
    return a

def solve(G, H, m):
    LG = nx.line_graph(G)
    LH = nx.line_graph(H)
    GN, HN = are_isomorphic(G, H, linear=False)
    GE, HE = are_isomorphic(LG, LH, linear=True)

    gtn = topological_nodes(G, GN)
    htn = topological_nodes(H, HN)
    gte = topological_edges(G, GE)
    hte = topological_edges(H, HE)

    inv_map = {v: k for k, v in m.items()}

    return gtn, [inv_map[k] for k in htn], gte, [(inv_map[i], inv_map[j]) for i, j in hte]


def assume(gtn, gte, htn, hte):
    assumptions = defaultdict(list)
    shift = 0
    for hn in htn:
        assumptions[hn] = gtn[shift:len(gtn) - len(htn) + shift + 1]
        shift += 1

    essumptions = defaultdict(list)
    shift = 0
    for he in hte:
        essumptions[he] = gte[shift:len(gte) - len(hte) + shift + 1]
        shift += 1

    try:
     while True:
        changed = False
        for e in essumptions.keys():
            for j in essumptions[e]:
                a, b = e
                x, y = j
                if x not in assumptions[a] or y not in assumptions[b]:
                    essumptions[e].remove(j)
                    changed = True

        for e in essumptions.keys():
            a, b = zip(*essumptions[e])
            a = set(a)
            b = set(b)
            x, y = e
            fg = set(assumptions[x]) - a
            hg = set(assumptions[y]) - b

            if len(fg) > 0:
                for f in fg:
                    assumptions[x].remove(f)
                    changed = True

            if len(hg) > 0:
                for h in hg:
                    assumptions[y].remove(h)
                    changed = True

        if not changed:
            break

    except:
        print()

    print(assumptions)
    print(essumptions)


if __name__ == '__main__':
    N = 30
    D = 5

    G, H, m = generate_problem_sub(D, N, isomorphic=True)
    gtn, htn, gte, hte = solve(G, H, m)

    print(gtn)
    print(htn)
    print(gte)
    print(hte)

    assume(gtn, gte, htn, hte)
