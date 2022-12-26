import random
import time
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

import itertools
import leidenalg
import igraph as ig
import networkx as nx

from collections import OrderedDict, Counter

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
    B = A.copy()

    rm = random.sample(list(B.nodes()), 10)

    k = list(B.nodes)
    random.shuffle(k)
    H = nx.Graph()
    H.add_nodes_from(k)
    H.add_edges_from(B.edges(data=True))
    B = H

    B.remove_nodes_from(rm)
    A.remove_nodes_from(list(nx.isolates(A)))
    B.remove_nodes_from(list(nx.isolates(B)))

    if not isomorphic:
        B = nx.double_edge_swap(B)
    if isomorphic != nx.isomorphism.GraphMatcher(A, B).subgraph_is_isomorphic():
        return generate_problem_sub(degree, size, isomorphic)

    return A, B


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
    e = Counter.most_common(VN)

    l = {}
    for x,y in e:
      a,b = x
      if a not in l:
          l[a] = (tuple([b]), y)
      else:
          l[a] = (tuple(list(l[a][0]) + [b]), l[a][1] + y)

    s = {}
    for x,y in l.items():
        s[(x, y[0])] = y[1]

    a,b = zip(*Counter.most_common(s))
    return a

def solve(G, H):
    LG = nx.line_graph(G)
    LH = nx.line_graph(H)
    GN, HN = are_isomorphic(G, H, linear=False)
    GE, HE = are_isomorphic(LG, LH, linear=True)

    gtn = topological_nodes(G, GN)
    htn = topological_nodes(H, HN)
    gte = topological_edges(G, GE)
    hte = topological_edges(H, HE)

    return gtn, htn, gte, hte


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

    print(assumptions)
    print(essumptions)


def generate_problem(degree, size, isomorphic):
    A = nx.random_regular_graph(degree, size)
    #A = random_k_out_graph(size, 3 , 0.6)
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
        #B = directed_edge_swap(B)
    if isomorphic != nx.is_isomorphic(A, B):
        return generate_problem(degree, size, isomorphic)

    return A, B

if __name__ == '__main__':
    N = 30
    D = 5

    G, H = generate_problem_sub(D, N, isomorphic=True)
    gtn, htn, gte, hte = solve(G, H)

    print(gtn)
    print(htn)
    print(gte)
    print(hte)

    #assume(gtn, gte, htn, hte)
