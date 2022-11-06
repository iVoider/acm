from collections import Counter

import numpy as np


def sandpile(g, sinks, steps = 1000):
    nodes = list(g.nodes)
    sand = dict.fromkeys(nodes, len(nodes))

    for sink in sinks:
        sand[sink] = 0

    waste = dict.fromkeys(sinks, 0)

    for s in sinks:
        sand[s] = 0

    for step in range(0, steps):
        node = max(sand, key=sand.get)
        if sand[node] > 0:
                to = -1
                to_value = len(nodes) ** 10
                for d in g.neighbors(node):
                    if sand[d] < to_value:
                        to = d
                        to_value = sand[d]
                if to != -1:
                    if to in sinks:
                        waste[to] += 1
                    else:
                     sand[to] += 1
                    sand[node] -= 1

    return (np.min(list(sand.values())) + 1) / (np.max(list(sand.values())) + 1)



