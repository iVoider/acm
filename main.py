import itertools
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

import random
import leidenalg
import igraph as ig
import networkx as nx
from collections import OrderedDict
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
        val = 1.0 - leidenalg.find_partition(g, partition_type=leidenalg.ModularityVertexPartition, seed=seed).modularity
        result[s] = val
        return val

    return f


def are_isomorphic(G, H, linear):
    def optimization(g, linear,seed):
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


def generate_problem(degree, size, isomorphic=True):
    A = nx.random_regular_graph(degree, size)
    node_mapping = dict(zip(A.nodes(), sorted(A.nodes(), key=lambda k: random.random())))
    B = nx.relabel_nodes(A, node_mapping)

    if not isomorphic:
        B = nx.double_edge_swap(B)

    return A, B, node_mapping


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

    B.remove_nodes_from(random.sample(list(B.nodes()), 10))
    A.remove_nodes_from(list(nx.isolates(A)))
    B.remove_nodes_from(list(nx.isolates(B)))

    if not isomorphic:
        B = nx.double_edge_swap(B)
    if isomorphic != nx.isomorphism.GraphMatcher(A, B).subgraph_is_isomorphic():
        return generate_problem_sub(degree, size, isomorphic)

    return A, B, node_mapping


if __name__ == '__main__':
     N = 30
     D = 5

     G, H, m = generate_problem_sub(D, N, isomorphic=False)
     LG = nx.line_graph(G)
     LH = nx.line_graph(H)
     GN, HN = are_isomorphic(G, H, linear=False)
     GE, HE = are_isomorphic(LG, LH, linear=True)

     m = {v: k for k, v in m.items()}

     VN = dict.fromkeys(G.nodes, (0, 0))
     for k, v in GN.items():
        for e in k:
            VN[e] = (VN[e][0] + v, VN[e][1] + 1)
     for k in VN:
        VN[k] = VN[k][0] / VN[k][1]
     print(list(VN.keys()))


     VN = dict.fromkeys(H.nodes, (0, 0))
     for k, v in HN.items():
         for e in k:
             VN[e] = (VN[e][0] + v, VN[e][1] + 1)
     for k in VN:
         VN[k] = VN[k][0] / VN[k][1]

     print([m[i] for i in VN.keys()])

     VN = dict.fromkeys(G.edges, (0, 0))
     for k, v in GE.items():
         for e in k:
             VN[e] = (VN[e][0] + v, VN[e][1] + 1)
     for k in VN:
       VN[k] = VN[k][0] / VN[k][1]
     A = list(VN.keys())
     print(A)

     VN = dict.fromkeys(H.edges, (0, 0))
     for k, v in HE.items():
         for e in k:
             VN[e] = (VN[e][0] + v, VN[e][1] + 1)
     for k in VN:
         VN[k] = VN[k][0] / VN[k][1]

     B = [(m[i[0]], m[i[1]]) for i in VN.keys()]
     print(B)
