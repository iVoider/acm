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
    while len(cur) > M:
      mx = dict.fromkeys(cur - choice, 0)
      for x in cur - choice:
            mx[x] = (sum([ap[i] for i in G.neighbors(x)]) + ap[x])
      m = max(mx, key=mx.get)
      choice.add(m)
      inc = set(G.neighbors(m))
      inc.add(m)
      cur &= inc
      todel = set()
      for c in set(G.vs.indices) - cur:
          todel |= set(G.vs[c].incident())
      G.delete_edges(todel)
      ap = aprob(G)

    print(G.subgraph(cur).clique_number(), M)
    break
