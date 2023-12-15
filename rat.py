def aprob(Gr):
  stationary_distrib = list()
  for v in Gr.vs:
    stationary_distrib.append(Gr.degree(v) / (2 * G.ecount()))
  return stationary_distrib

def bry(G):
  best = dict()
  H = G.copy()
  apn = aprob(G)

  for v in G.vs.indices:
      best[v] = (sum([apn[i] for i in G.neighbors(v)]) + apn[v])

  for v in sorted(best, key = best.get, reverse=True):
    ap = apn.copy()
    G = H.copy()
    cur = set(G.neighbors(v))
    cur.add(v)
    choice = set({v})
    for zi in range(1, M + 1):
      mx = dict.fromkeys(cur, 0)
      for x in cur:
            mx[x] = (sum([ap[i] for i in G.neighbors(x)]) + ap[x])
      choice = set(sorted(mx, key=mx.get, reverse=True)[0:zi])
      inc = cur.copy()
      for c in choice:
        p = set(G.neighbors(c))
        p.add(c)
        inc &= p

      cur = inc.copy()

      todel = set()
      for x in set(G.vs.indices) - cur:
          todel |= set(G.incident(x))

      G.delete_edges(todel)
      ap = aprob(G)

    print(G.subgraph(choice).clique_number(), M)
    break
