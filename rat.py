def aprob(Gr):
  stationary_distrib = list()
  for v in Gr.vs:
    stationary_distrib.append(Gr.degree(v) / (2 * G.ecount()))
  return stationary_distrib

def recalc(G):
  ap = aprob(G)

  res = dict([(v["_nx_name"], list()) for v in G.vs])

  all_sat = dict()

  for v in G.vs:
    nv = [G.vs[n]["_nx_name"] for n in G.neighbors(v)]
    nv.append(v["_nx_name"])
    all_sat[v["_nx_name"]] = set(nv)

  prios = dict()
  for v in G.vs.indices:
    global_vals = dict()

    z = sum([ap[x] for x in G.neighbors(v)])
    for n in G.neighbors(v):
      global_vals[G.vs[n]["_nx_name"]] = ap[n] / z
    res[G.vs[v]["_nx_name"]] = global_vals
    prios[v] = z + ap[v].real
  return sorted(prios, key=prios.get, reverse = True), all_sat, res

def sald(G):
  orig = dict([(v["_nx_name"], v.index) for v in G.vs])

  p, all_sat0, res0 = recalc(G)
  for vp in p:
    all_sat, res = all_sat0, res0
    v = G.vs[vp]["_nx_name"]
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
      # _, all_sat, res = recalc(G.subgraph([orig[c] for c in cur]))
    zoo = G.subgraph([orig[c] for c in cur]).clique_number()
    print(zoo)
    if zoo == M:
      print("Booba")
      break
