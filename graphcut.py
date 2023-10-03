N = 10
M = N * 4

for _ in range(0, 10):
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

  found = False
  for v in G.vs.indices:
    sub = all_sat[v]
    cur = set({v})
    while len(sub) > M:
     ls = sub - cur
     val = dict.fromkeys(ls, 0)
     for n in ls:
      for z in sub - set({n}):
         if n in res[z]:
            val[n] += res[z][n]
     mx = max(val, key=val.get)
     cur.add(mx)
     sub = all_sat[v]
     for c in cur:
      sub &= all_sat[c]
    
    
    if len(sub) == M:
      h = G.subgraph(sub)
      if h.density() == 1.0:
        found = True
        break

  if not found:
     print('Shit')
