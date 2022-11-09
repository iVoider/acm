import random

import networkx as nx
from networkx import random_regular_graph, double_edge_swap, gnp_random_graph, random_k_out_graph
from networkx.utils import py_random_state


@py_random_state(3)
def directed_edge_swap(G, *, nswap=1, max_tries=100, seed=None):
    """Swap three edges in a directed graph while keeping the node degrees fixed.
    A directed edge swap swaps three edges such that a -> b -> c -> d becomes
    a -> c -> b -> d. This pattern of swapping allows all possible states with the
    same in- and out-degree distribution in a directed graph to be reached.
    If the swap would create parallel edges (e.g. if a -> c already existed in the
    previous example), another attempt is made to find a suitable trio of edges.
    Parameters
    ----------
    G : DiGraph
       A directed graph
    nswap : integer (optional, default=1)
       Number of three-edge (directed) swaps to perform
    max_tries : integer (optional, default=100)
       Maximum number of attempts to swap edges
    seed : integer, random_state, or None (default)
        Indicator of random number generation state.
        See :ref:`Randomness<randomness>`.
    Returns
    -------
    G : DiGraph
       The graph after the edges are swapped.
    Raises
    ------
    NetworkXError
        If `G` is not directed, or
        If nswap > max_tries, or
        If there are fewer than 4 nodes in `G`
    NetworkXAlgorithmError
        If the number of swap attempts exceeds `max_tries` before `nswap` swaps are made
    Notes
    -----
    Does not enforce any connectivity constraints.
    The graph G is modified in place.
    References
    ----------
    .. [1] Erdős, Péter L., et al. “A Simple Havel-Hakimi Type Algorithm to Realize
           Graphical Degree Sequences of Directed Graphs.” ArXiv:0905.4913 [Math],
           Jan. 2010. https://doi.org/10.48550/arXiv.0905.4913.
           Published  2010 in Elec. J. Combinatorics (17(1)). R66.
           http://www.combinatorics.org/Volume_17/PDF/v17i1r66.pdf
    .. [2] “Combinatorics - Reaching All Possible Simple Directed Graphs with a given
           Degree Sequence with 2-Edge Swaps.” Mathematics Stack Exchange,
           https://math.stackexchange.com/questions/22272/. Accessed 30 May 2022.
    """
    if nswap > max_tries:
        raise nx.NetworkXError("Number of swaps > number of tries allowed.")
    if len(G) < 4:
        raise nx.NetworkXError("Graph has less than four nodes.")

    # Instead of choosing uniformly at random from a generated edge list,
    # this algorithm chooses nonuniformly from the set of nodes with
    # probability weighted by degree.
    tries = 0
    swapcount = 0
    keys, degrees = zip(*G.degree())  # keys, degree
    cdf = nx.utils.cumulative_distribution(degrees)  # cdf of degree
    discrete_sequence = nx.utils.discrete_sequence

    while swapcount < nswap:
        # choose source node index from discrete distribution
        start_index = discrete_sequence(1, cdistribution=cdf, seed=seed)[0]
        start = keys[start_index]
        tries += 1

        if tries > max_tries:
            msg = f"Maximum number of swap attempts ({tries}) exceeded before desired swaps achieved ({nswap})."
            raise nx.NetworkXAlgorithmError(msg)

        # If the given node doesn't have any out edges, then there isn't anything to swap
        if G.out_degree(start) == 0:
            continue
        second = seed.choice(list(G.succ[start]))
        if start == second:
            continue

        if G.out_degree(second) == 0:
            continue
        third = seed.choice(list(G.succ[second]))
        if second == third:
            continue

        if G.out_degree(third) == 0:
            continue
        fourth = seed.choice(list(G.succ[third]))
        if third == fourth:
            continue

        if (
            third not in G.succ[start]
            and fourth not in G.succ[second]
            and second not in G.succ[third]
        ):
            # Swap nodes
            G.add_edge(start, third)
            G.add_edge(third, second)
            G.add_edge(second, fourth)
            G.remove_edge(start, second)
            G.remove_edge(second, third)
            G.remove_edge(third, fourth)
            swapcount += 1

    return G

def generate_problem(degree, size, isomorphic, directed=True):
    A = random_regular_graph(degree, size)
    #A = random_k_out_graph(size, 3 , 0.6)
    node_mapping = dict(zip(A.nodes(), sorted(A.nodes(), key=lambda k: random.random())))
    B = nx.relabel_nodes(A, node_mapping)

    A.remove_nodes_from(list(nx.isolates(A)))
    B.remove_nodes_from(list(nx.isolates(B)))

    if not isomorphic:
        B = double_edge_swap(B)
        #B = directed_edge_swap(B)
    if isomorphic != nx.is_isomorphic(A, B):
        return generate_problem(degree, size, isomorphic)

    return A, B, node_mapping

