import difflib
import itertools
import math
import warnings

import numpy as np

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

def test_arbitrage(s, d, add = True):
    val = 0
    for p in itertools.permutations(s):
        p = list(p)
        if add:
         p.append(p[0])

        v = 1
        for x, y in zip(p, p[1:]):
            v *= d[x][y]

        val = max(val, v)
    return val


def f_wrapper(d, result):

    def f(s):
        s = tuple(sorted(s))
        if len(s) < 2:
            return 0

        if s in result:
            return result[s]

        val = test_arbitrage(s, d)

        result[s] = val
        return val

    return f


def solve(d):
    result = dict()
    queyranne(f_wrapper(d, result), list(range(0, len(d))))
    return result


if __name__ == '__main__':

 d = [[1,1.35,0.93,0.83,7.85,0.92,134.1,1.45,82.82,6.87],
      [0.74,1,0.69,0.62,5.82,0.69,99.53,1.08,61.46,5.1],
      [1.07,1.44,1,0.89,8.41,0.99,143.71,1.56,88.75,7.36],
      [1.2,1.62,1.12,1,9.45,1.11,161.52,1.75,99.74,8.27],
      [0.13,0.17,0.12,0.11,1,0.12,17.09,0.19,10.55,0.88],
      [1.08,1.46,1.01,0.9,8.48,1,144.99,1.57,89.54,7.43],
      [0.01,0.01,0.01,0.01,0.06,0.01,1,0.01,0.62,0.05],
      [0.69,0.93,0.64,0.57,5.4,0.64,92.25,1,56.97,4.72],
      [0.01,0.02,0.01,0.01,0.09,0.01,1.62,0.02,1,0.08],
      [0.15,0.2,0.14,0.12,1.14,0.13,19.53,0.21,12.06,1]
      ]

 k = solve(d)
 print(Counter.most_common(k))

print(test_arbitrage([8,7,8,6,0], d, False))
