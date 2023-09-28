N = 5
M = N * 4

for _ in range(0, 100):
  f = gen_sat(N, M, random_cnf(N))
  g, _ = sat_to_clique(f)
  G = ig.Graph.from_networkx(g)

  ap = aprob(G)

  res = dict([(v, list()) for v in G.vs.indices])

  for v in G.vs.indices:
    local = set(G.neighbors(v) + [v])
    global_list = list()
    for _ in range(0, M ** 2):
     cur = v
     while True:
      rng = np.random.default_rng()
      probs = [ap[x].real for x in G.neighbors(cur)]
      probs.append(1.0 - sum(probs))
      next = np.random.choice(G.neighbors(cur) + [cur], p=probs)
      if next in local:
        cur = next
        global_list.append(cur)
      else:
        break
    r = dict(Counter(global_list).most_common())
    res[v] = r
  
  found = False
  for v in G.vs.indices:
    sub = set(G.neighbors(v) + [v])
    cur = set({v})
    while len(sub) > M:
     val = dict.fromkeys(sub - cur, 0)
     for n in sub - cur:
      for z in sub - set({n}):
         if n in res[z]:
            val[n] += res[z][n]
     mx = max(val, key=val.get)
     cur.add(mx)
     sub = set(G.neighbors(v) + [v])
     for c in cur:
      sub &= set(G.neighbors(c) + [c])
    if G.subgraph(sub).density() == 1.0 and G.subgraph(sub).vcount() == M:
      found = True
      break
  print(found)
