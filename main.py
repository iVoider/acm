def lnaprob(Gr):
    stationary_distrib = [0.0] * Gr.ecount()
    enum = 0
    for e in Gr.es:
           r = len(Gr.neighbors(e.source)) + len(Gr.neighbors(e.target)) - 2
           stationary_distrib[e.index] = r
           enum += r

    for i,e in enumerate(stationary_distrib):
        stationary_distrib[i] = e / ( enum)

    return stationary_distrib

def stat(transition):
    transition_matrix = np.array(transition)
    transition_matrix = transition_matrix/transition_matrix.sum(axis=1,keepdims=1)
    transition_matrix_transp = transition_matrix.T
    eigenvals, eigenvects = np.linalg.eig(transition_matrix_transp)
    close_to_1_idx = np.isclose(eigenvals, 1)
    target_eigenvect = eigenvects[:, close_to_1_idx]
    target_eigenvect = target_eigenvect[:, 0]
    stationary_distrib = target_eigenvect / sum(target_eigenvect)
    return [a.real for a in stationary_distrib]

def create_t(Gr, vals):
    transition = []
    for v in Gr.vs:
        cur = [0] * Gr.vcount()
        for n in v.neighbors():
            cur[n.index] = (vals[Gr.get_eid(v.index,n)])
        transition.append(cur)
    return transition

def aprob(Gr):
  stationary_distrib = list()
  for v in Gr.vs:
    stationary_distrib.append(Gr.degree(v) / (2 * Gr.ecount()))
  return stationary_distrib
