import itertools
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

import random
import leidenalg
import igraph as ig
import networkx as nx
from collections import OrderedDict, deque, defaultdict
import time


def queyranne(F, V):
    def Fnew(a):
        r = []
        for x in a:
            r += S[x - 1]
        return F(r)

    n = len(V)
    S = [[x] for x in V]
    inew = OrderedDict()
    for x in range(1, n + 1):
        inew[x] = x
    for h in range(n - 1):
        pendentpair(Fnew, inew)


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
        val = 1.0 - ig.community._community_leiden(g).modularity
        result[s] = val
        return val

    return f


def are_isomorphic(G, H, linear):
    def optimization(g, linear, seed):
        result = {}
        G = ig.Graph.from_networkx(g)
        queyranne(f_wrapper(G, result, seed), list(range(0, g.number_of_nodes())))
        ret = {}
        for key, value in result.items():
            ret[tuple(sorted(
                [(int(G.vs[i]["_nx_name"]) if not linear else tuple(G.vs[i]["_nx_name"])) for i in key]))] = value
        return ret

    seed = random.getrandbits(16)
    A = optimization(G, linear, seed)
    B = optimization(H, linear, seed)
    return A, B


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
    A.remove_nodes_from(list(nx.isolates(A)))
    B.remove_nodes_from(list(nx.isolates(B)))

    if not isomorphic:
        B = nx.double_edge_swap(B)
    if isomorphic != nx.isomorphism.GraphMatcher(A, B).subgraph_is_isomorphic():
        return generate_problem_sub(degree, size, isomorphic)

    return A, B, node_mapping, rm


def topological_nodes(G, GN):
    VN = dict.fromkeys(G.nodes, (0, 0))
    for k, v in GN.items():
        for e in k:
            VN[e] = (VN[e][0] + v, VN[e][1] + 1)
    for k in VN:
        VN[k] = VN[k][0] / VN[k][1]
    return list(VN.keys())


def topological_edges(G, GE):
    VN = dict.fromkeys(G.edges, (0, 0))
    for k, v in GE.items():
        for e in k:
            VN[e] = (VN[e][0] + v, VN[e][1] + 1)
    for k in VN:
        VN[k] = VN[k][0] / VN[k][1]
    return list(VN.keys())


def traverse(L_N, L_E, G):
    Q = deque()
    used = dict.fromkeys(L_N, False)
    S = list()

    res = defaultdict(list)
    for i, j in L_E:
        res[i].append(j)

    for i, j in L_E:
        res[j].append(i)

    for k in res:
        v = list()
        for p in L_E:
            for t in res[k]:
                if (t, k) == p or (k, t) == p:
                    v.append(t)
        res[k] = v

    Q.append(L_E[0])
    used[L_E[0]] = 1
    label = 2
    N = list()
    J = list()
    while Q:
        vi, vj = Q.popleft()
        for v in (res[vi][res[vi].index(vj):] + res[vi][:res[vi].index(vj)]):
            if not used[v]:
                Q.append((v, vi))
                used[v] = label
                label += 1
            N.append(v)
            J.append(leidenalg.find_partition(ig.Graph.from_networkx(G.subgraph(N)),
                                              leidenalg.ModularityVertexPartition).modularity)
            S.append(used[v])

    return S, N, J


def traverse_2(L_N, L_E, G, CJ):
    Q = deque()
    used = dict.fromkeys(L_N, False)
    S = list()

    res = defaultdict(list)
    for i, j in L_E:
        res[i].append(j)

    for i, j in L_E:
        res[j].append(i)

    for k in res:
        v = list()
        for p in L_E:
            for t in res[k]:
                if (t, k) == p or (k, t) == p:
                    v.append(t)
        res[k] = v

    Q.append(L_E[0])
    used[L_E[0]] = 1
    label = 2
    N = list()
    J = list()
    cur = 0
    while Q:
        vi, vj = Q.popleft()
        for v in (res[vi][res[vi].index(vj):] + res[vi][:res[vi].index(vj)]):
            KN = N + [v]
            r = leidenalg.find_partition(ig.Graph.from_networkx(G.subgraph(KN)), leidenalg.ModularityVertexPartition).modularity
            if str(r) == "nan" and str(CJ[cur]) == "nan":
                chk = True
            else:
             chk = r == CJ[cur]

            if not used[v] and chk:
                Q.append((v, vi))
                used[v] = label
                label += 1
            if chk:
             N.append(v)
             J.append(leidenalg.find_partition(ig.Graph.from_networkx(G.subgraph(N)),
                                              leidenalg.ModularityVertexPartition).modularity)
             S.append(used[v])
             cur += 1
             if len(J) == len(CJ):
                 return S,N,J

    return S, N, J


def solve(G, H, mi):
    LG = nx.line_graph(G)
    LH = nx.line_graph(H)
    GN, HN = are_isomorphic(G, H, linear=False)
    GE, HE = are_isomorphic(LG, LH, linear=True)

    m = {v: k for k, v in mi.items()}

    gtn = topological_nodes(G, GN)
    htn = topological_nodes(H, HN)
    gte = topological_edges(G, GE)
    hte = topological_edges(H, HE)

    print(gtn)
    print([m[i] for i in htn])
    print(gte)
    print([(m[i[0]], m[i[1]]) for i in hte])

    print()

    C1, N1, J1 = traverse(htn, hte, H)
    C2, N2, J2 = traverse_2(gtn, gte, G, J1)
    print(C1)
    print(C2)
    print()
    print([m[i] for i in N1])
    print(N2)
    print()
    print(J1)
    print(J2)


if __name__ == '__main__':
    N = 25
    D = 4

    G, H, mi, rm = generate_problem_sub(D, N, isomorphic=True)
    solve(G, H, mi)
