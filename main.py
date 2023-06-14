def diffactor(G):

    zol = dict.fromkeys(G.vs.indices, 0)

    for e in G.es:

      x = e.source
      y = e.target

      teamX = set(G.neighbors(x)) - set(G.neighbors(y))
      teamY = set(G.neighbors(y)) - set(G.neighbors(x))
      teamM = set(G.neighbors(x)) & set(G.neighbors(y))
      teamX.add(x)
      teamY.add(y)
      if x in teamY:
       teamY.remove(x)

      if y in teamX:
       teamX.remove(y)

      if y in teamM:
       teamM.remove(y)

      if x in teamM:
       teamM.remove(x)

      prev = len(teamX) + len(teamY) + len(teamM)
      total = prev

      while prev < G.vcount():
       tx = teamX.copy()
       ty = teamY.copy()
       for peace in teamM.copy():
        n = set(G.neighbors(peace))
        a = len(n & tx)
        b = len(n & ty)
        if  a > b:
          teamX.add(peace)
          teamM.remove(peace)
        elif a < b:
          teamY.add(peace)
          teamM.remove(peace)

       tx = teamX.copy()
       ty = teamY.copy()
       for free in set(G.vs.indices) - (tx | ty | teamM):
        n = set(G.neighbors(free))
        a = len(n & tx)
        b = len(n & ty)
        if a > b:
          teamX.add(free)
        elif a < b:
          teamY.add(free)

       total = len(teamX) + len(teamY) + len(teamM)
       if total == prev:
        break
       else:
        prev = total

      if len(teamX) > len(teamY):
        zol[x] += 1
      elif len(teamX) < len(teamY):
        zol[y] += 1

    return zol


def f_wrapper(G, size, result):
    def f(V, s, params=None):

        if type(s) == int:
            s = [list(G.vs.indices)[s]]
        else:
            s = list(itertools.chain(*(i if isinstance(i, list) else (i,) for i in s)))
            s = [list(G.vs.indices)[n] for n in s]

        s = tuple(sorted(s))

        if s in result:
            return result[s]

        p = set(G.vs.indices) - set(s)
        if len(p) == 0:
          result[s] = -100
        else:
         result[s] = (sum(diffactor(G.subgraph(p)).values()) / len(p)) * -1
        return result[s]

    return f


def solve(G, size):
    result = {}

    transition = []

    for v in G.vs:
        nxt = [0.0] * G.vcount()
        for n in v.neighbors():
            nxt[n.index] = 1.0
        transition.append(nxt)

    Q = np.array(transition)
    io = QUEYRANNE(Q, f_wrapper(G, size, result))
    return io
