N = 300
M = N * 4

for _ in range(0, 1):
  f = gen_sat(N, M, random_cnf(N))
  g, _ = sat_to_clique(f)
  G = ig.Graph.from_networkx(g)

  ap = aprob(G)

  res = dict([(v, list()) for v in G.vs.indices])

  all_sat = dict()

  for v in G.vs.indices:
    nv = G.neighbors(v)
    nv.append(v)
    all_sat[v] = set(nv)

  for v in G.vs.indices:
    global_vals = dict()
    z = sum([ap[x].real for x in G.neighbors(v)])
    for n in G.neighbors(v):
      global_vals[n] = ap[n].real * 10 /  z * 10
    res[v] = global_vals

  mxs = set()
  for v in G.vs.indices:
    cur = all_sat[v]
    choice = set({v})
    while len(cur) > M:
      mx = dict.fromkeys(cur - choice, 0)
      for x, y in itertools.product(cur - choice, choice):
         if y in res:
            mx[x] += res[y][x]
      m = max(mx, key=mx.get)
      choice.add(m)
      cur &= all_sat[m]
    mxs.add(G.subgraph(cur).clique_number())
    print(max(mxs), M)
    break
