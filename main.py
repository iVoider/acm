def diffactor(G):

    zol = dict([(v, list()) for v in G.vs.indices])

    for x,y,z in itertools.combinations(G.vs.indices, 3):

      if G.get_eid(x,y,error=False) != -1 or G.get_eid(x,z,error=False) != -1 or G.get_eid(z,y, error=False) != -1:
        continue

      teamX = set(G.neighbors(x)) - set(G.neighbors(y)) - set(G.neighbors(z))
      teamY = set(G.neighbors(y)) - set(G.neighbors(x)) - set(G.neighbors(z))
      teamZ = set(G.neighbors(z)) - set(G.neighbors(y)) - set(G.neighbors(x))
      teamM = set(G.neighbors(x)) & set(G.neighbors(y)) & set(G.neighbors(z))

      teamX.add(x)
      teamY.add(y)
      teamZ.add(z)

      if x in teamY:
       teamY.remove(x)

      if y in teamX:
       teamX.remove(y)

      if y in teamZ:
       teamZ.remove(y)

      if x in teamZ:
       teamZ.remove(x)

      if z in teamX:
        teamX.remove(z)

      if z in teamY:
        teamY.remove(z)

      if y in teamM:
       teamM.remove(y)

      if z in teamM:
       teamM.remove(z)

      if x in teamM:
       teamM.remove(x)

      prev = len(teamX) + len(teamY) + len(teamZ) + len(teamM)
      total = prev

      while prev < G.vcount():
       tx = teamX.copy()
       ty = teamY.copy()
       tz = teamZ.copy()

       for peace in teamM.copy():
        n = set(G.neighbors(peace))
        a = len(n & tx)
        b = len(n & ty)
        c = len(n & tz)

        tr = [a, b, c]
        if len(set(tr)) == 3:
         if tr.index(max(a, b, c)) == 0:
          teamX.add(peace)
          teamM.remove(peace)
         elif tr.index(max(a, b, c)) == 1:
          teamY.add(peace)
          teamM.remove(peace)
         elif tr.index(max(a, b, c)) == 2:
          teamZ.add(peace)
          teamM.remove(peace)

       tx = teamX.copy()
       ty = teamY.copy()
       tz = teamZ.copy()

       for free in set(G.vs.indices) - (tx | ty | tz | teamM):
        n = set(G.neighbors(free))
        a = len(n & tx)
        b = len(n & ty)
        c = len(n & tz)

        tr = [a, b, c]
        if len(set(tr)) == 3:
         if tr.index(max(a, b, c)) == 0:
          teamX.add(free)
         elif tr.index(max(a, b, c)) == 1:
          teamY.add(free)
         elif tr.index(max(a, b, c)) == 2:
          teamZ.add(free)

       total = len(teamX) + len(teamY) + len(teamM) + len(teamZ)
       if total == prev:
        break
       else:
        prev = total

      zol[x] = len(teamX) - len(teamM)
      zol[y] = len(teamY) - len(teamM)
      zol[z] = len(teamZ) - len(teamM)

    cur = dict.fromkeys([-1,-2,-3,-4,-5,1,2,3,4,5], 0)

    for v in G.vs.indices:
      frcr = 0
      frcr += zol[v]
      for n in G.neighbors(v):
        frcr += zol[n]
      cur[G.vs[v]["_nx_name"][0]] += frcr

    return np.product([abs(cur[1]) - abs(cur[-1]), abs(cur[2]) - abs(cur[-2]), abs(cur[3]) - abs(cur[-3]), abs(cur[4]) - abs(cur[-4]), abs(cur[5]) - abs(cur[-5])])

yes = set()
no = set()

M = 20
N = 5

for i in range(0, 1000):
 f = gen_sat(N, M, random_cnf(N))
 g, _ = sat_to_clique(f)
 G = ig.Graph.from_networkx(g)
 d = diffactor(G)
 yes.add(d)

for i in range(0, 1000):
 f = gen_unsat(N, M)
 g, _ = sat_to_clique(f)
 G = ig.Graph.from_networkx(g)
 d = diffactor(G)
 no.add(d)

print(len(yes), len(no))

print(len(yes & no))

print(len(yes - no),len(no - yes) )

dx = sorted(yes) + sorted(no)
dy = [1] * len(yes) + [0] * len(no)
dc = [1] * len(yes) + [0] * len(no)
import plotly.express as px
fig = px.scatter(x=dx, y=dy,color=dc)
fig.show()
