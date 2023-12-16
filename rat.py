def aprob(Gr):
  stationary_distrib = list()
  for v in Gr.vs:
    stationary_distrib.append(Gr.degree(v) / (2 * Gr.ecount()))
  return stationary_distrib

def bry(G):
  best = dict()
  H = G.copy()
  apn = aprob(G)

  fnd = False

  for v in G.vs.indices:
      best[v] = (apn[v], sum([apn[i] for i in G.neighbors(v)]))

  for v in sorted(best, key = best.get, reverse=True):
    ap = apn.copy()
    G = H.copy()
    cur = set(G.neighbors(v))
    cur.add(v)
    choice = set({v})
    for zi in range(1, M + 1):
      mx = dict.fromkeys(cur, 0)
      for x in cur:
            mx[x] = (ap[x], sum([ap[i] for i in G.neighbors(x)]))
      choice = set(sorted(mx, key=mx.get, reverse=True)[0:zi])
      inc = cur.copy()
      for c in choice:
        p = set(G.neighbors(c))
        p.add(c)
        inc &= p
        if len(inc) < len(cur):
           break

      cur = inc.copy()

      todel = set()
      for x in set(G.vs.indices) - cur:
          todel |= set(G.incident(x))

      G.delete_edges(todel)
      while True:
        changed = False
        for v in G.vs.indices:
          if G.degree(v) < (M - 1) and G.degree(v) > 0:
             G.delete_edges(G.incident(v))
             changed = True
        if not changed:
           break

      ap = aprob(G)
    t = G.subgraph(choice).clique_number()
    if t == M:
       fnd = True
       break

  assert(fnd)
