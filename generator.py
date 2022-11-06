import itertools
import random
import networkx as nx

def random_edge(nodes, ignore, used_edges, start):
    while True:
        choice = random.choice(nodes)
        if choice != ignore and (start, choice) not in used_edges:
            return start, choice


def generate_rand_digraphs(N, sparse=0.3, dropout=8, isomorphic=True):
    A = nx.MultiDiGraph()
    B = nx.MultiDiGraph()
    nodesA = list(range(0, N))
    nodesB = nodesA.copy()
    random.shuffle(nodesB)
    mapping = dict(zip(nodesA, nodesB))
    ins = dict(zip(nodesB, [set() for _ in range(N)]))
    outs = dict(zip(nodesB, [set() for _ in range(N)]))

    insert = []

    while len(A.nodes) < N:
        x = random.randint(0, N - 1)
        y = random.randint(0, N - 1)

        if random.random() < sparse:
            A.add_edge(x, y)
            insert.append((mapping[x], mapping[y]))
            ins[mapping[x]].add(mapping[y])
            outs[mapping[y]].add(mapping[x])

        if random.random() < sparse:
            A.add_edge(y, x)
            insert.append((mapping[y], mapping[x]))
            ins[mapping[y]].add(mapping[x])
            outs[mapping[x]].add(mapping[y])

    random.shuffle(insert)

    for edge in insert:
        B.add_edge(edge[0], edge[1])

    if not isomorphic:
        for x, y in itertools.combinations(nodesB, 2):
            if len(ins[x]) == len(ins[y]) and len(outs[x]) == len(outs[y]):
                if len(ins[x]) > 0 and len(ins[y]) > 0:
                    a = ins[x].pop()
                    b = ins[y].pop()
                    ins[x].add(b)
                    ins[y].add(a)
                    B.remove_edge(x, a)
                    B.add_edge(x, b)
                    B.remove_edge(y, b)
                    B.add_edge(y, a)
                    dropout -= 1
                elif len(outs[x]) > 0 and len(outs[y]) > 0:
                    a = outs[x].pop()
                    b = outs[y].pop()
                    outs[x].add(b)
                    outs[y].add(a)
                    B.remove_edge(a, x)
                    B.add_edge(b, x)
                    B.remove_edge(b, y)
                    B.add_edge(a, y)
                    dropout -= 1
            if dropout == 0:
                break

        if dropout != 0:
            return generate_rand_digraphs(N, sparse, dropout, isomorphic)

    return A, B, mapping