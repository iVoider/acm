

import random
from collections import OrderedDict, Counter

import numpy as np

from sat import gen_sat, random_cnf, gen_unsat, sat


def derandomizedSolve(f, n, m):
    # Assigning
    assignment = []
    clauseStatus = []

    for i in range(m):
        clauseStatus.append(-1)

    for i in range(n):

        trueCount = 0
        falseCount = 0

        for idx, clause in enumerate(f):

            if (clauseStatus[idx] == -1):
                if (i in clause):
                    trueCount += 1
                if (i + n in clause):
                    falseCount += 1

        assignment.append(trueCount > falseCount)
        for idx, clause in enumerate(f):
            if (i + (trueCount <= falseCount) * n in clause):
                clauseStatus[idx] = 1

    for i in range(n):
        assignment.append(not assignment[i])

    # Checking
    satisfiedClauses = []

    for i in range(m):
        if (assignment[f[i][0]] or assignment[f[i][1]] or assignment[f[i][2]]):
            satisfiedClauses.append(i)

    if len(satisfiedClauses) == 0:
        return 0

    return len(satisfiedClauses)


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

def cut_formula(f, s):
    fo = [f[i] for i in s]

    vars = {}
    vars_count = 1
    nf = []
    for c in fo:
        nc = []
        for v in c:
            if abs(v) not in vars:
                vars[abs(v)] = vars_count
                vars_count += 1
            nc.append(np.sign(v) * vars[abs(v)])
        nf.append(nc)
    return nf, vars_count, len(nf)


def f_wrapper(d, n, result):

    r = derandomizedSolve(d, n, len(d))

    def f(s):

        s = tuple(sorted(s))

        if s in result:
            return result[s]

        # cf, n , m = cut_formula(d, s)
        # val1 = derandomizedSolve(cf, n, m)


        cut = tuple(sorted(set(range(len(d))) - set(s)))
        cf, n, m = cut_formula(d, cut)
        val2 = derandomizedSolve(cf, n, m)

        result[s] = val2

        return result[s]

    return f


def solve(d, n):
    result = dict()
    ret = []
    wor = queyranne(f_wrapper(d, n, result), list(range(0, len(d))))
    for e, v in wor:
        ret.append(d[e[0]])
    return ret


if __name__ == '__main__':

    # Number of clauses
    m = 4 * 5

    # Number of variables
    n = 5

    total = [0] * m

    for j in range(0,10000):
     formula = gen_unsat(n, m)

     i = 0
     for c in solve(list(formula), n):
        if not sat(formula - {c}):
            total[i] += 1
        i += 1


    print(total)

