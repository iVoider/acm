N = 30
M = N * 4

rng = np.random.default_rng()

for _ in range(0, 1):
  f = gen_sat(N, M, random_cnf(N))
  g, _ = sat_to_clique(f)
  G = ig.Graph.from_networkx(g)

  ap = aprob(G)

  res = dict([(v, list()) for v in G.vs.indices])

  all_probs = dict()
  all_cat = dict()

  for v in G.vs.indices:
    nv = G.neighbors(v)
    pv = [ap[x].real for x in nv]
    pv.append(1.0 - sum(pv))
    all_probs[v] = pv
    nv.append(v)
    all_cat[v] = nv
  
  # all probs has some logic error unlike previous version
  for v in G.vs.indices:
    local = set(all_cat[v])
    global_vals = dict.fromkeys(G.vs.indices, 0)
    for _ in range(0, N * M):
     cur = v
     steps = 0
     while steps < N * M:
      next = random.choices(all_cat[cur], k = 1, weights=all_probs[cur])[0]
      if next in local:
        cur = next
        global_vals[cur] += 1
      else:
        break
      steps += 1
    res[v] = global_vals

  found = False
  for v in G.vs.indices:
    sub = set(all_cat[v])
    cur = set({v})
    while len(sub) > M:
     val = dict.fromkeys(sub - cur, 0)
     for n in sub - cur:
      for z in sub - set({n}):
         if n in res[z]:
            val[n] += res[z][n]
     mx = max(val, key=val.get)
     cur.add(mx)
     sub = set(all_cat[v])
     for c in cur:
      sub &= set(all_cat[c])
    h = G.subgraph(sub)
    if h.vcount() == M and h.density() == 1.0:
      found = True
      break
  if not found:
     print('Shit')
