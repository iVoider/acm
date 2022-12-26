import random
import time
import warnings

from networkx.algorithms.approximation import max_clique

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


def queyranne(F,V):
    def Fnew(a):
        r=[]
        for x in a:
            r += S[x - 1]
        return F(r)
    n = len(V)
    S = [[x] for x in V]
    s = []
    A = []
    inew = OrderedDict()
    for x in range(1,n+1):
        inew[x] = x
    minimum=float("inf")
    position_of_min=0
    for h in range(n-1):
        #Find a pendant pair
        [t,u] = pendentpair(Fnew,inew)
        #This gives a candidate solution
        A.append(S[u - 1].copy())
        s.append(Fnew({u}))
        if s[-1] < minimum:
            minimum = s[-1]
            position_of_min = len(s) - 1
        S[t - 1] += S[u - 1]
        del inew[u]
        for x in range(len(S[u - 1])):
            S[u - 1][x] *= -1
    vals = dict(zip([tuple(a) for a in A],s))
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
        val = 1.0 - leidenalg.find_partition(g, partition_type=leidenalg.ModularityVertexPartition, seed=seed).modularity
        result[s] = val
        return val

    return f

def is_complete_graph(G):
    N = len(G) - 1
    return not any(n in nbrdict or len(nbrdict)!=N for n, nbrdict in G.adj.items())

def is_subclique(G,nodelist):
    H = G.subgraph(nodelist)
    n = len(nodelist)
    return H.size() == n*(n-1)/2

def solve(G, H, s, linear):
    def optimization(g, linear, seed):
        result = {}
        G = ig.Graph.from_networkx(g)
        q = queyranne(f_wrapper(G, result, seed), list(range(0, g.number_of_nodes())))
        ret = {}
        for key, value in result.items():
          if len(key) >= s:
           ret[tuple(sorted(
           [(int(G.vs[i]["_nx_name"]) if not linear else tuple(G.vs[i]["_nx_name"])) for i in key]))] = value
        return q

    seed = random.getrandbits(16)
    x = optimization(G, linear, seed)
    y = optimization(nx.line_graph(G), True, seed)

    print(x)
    print(y[:1], y[-1:])

    return False


def generate_problem_clique(degree, size, contains):
    while True:
        try:
            A = nx.random_regular_graph(degree, size)
        except:
            A = None
        if A != None:
            break

    M = ig.Graph.from_networkx(A)

    largest = M.clique_number()
    print(M.largest_cliques())

    if contains:
     B = nx.complete_graph(largest)
    else:
     B = nx.complete_graph(largest + 1)

    k = list(B.nodes)
    random.shuffle(k)
    H = nx.Graph()
    H.add_nodes_from(k)
    H.add_edges_from(B.edges(data=True))
    B = H

    return A, B, largest

if __name__ == '__main__':
    N = 30
    D = 10

    G, H, s = generate_problem_clique(D, N, contains=True)
    print(s)
    print(solve(G, H, s, linear=False))
