"""
author: [Tamim Hussein]
studentnumber: [2502572]
"""

import matplotlib.pyplot as plt
import json
import networkx as nx


def extract_kmers(json_data, k):
    data = json.loads(json_data)
    kmers = list()
    for sequence in data.values():
        while len(sequence) >= k:
            kmers.append(sequence[:k])
            sequence = sequence[1:]
    return kmers


def construct_graph(json_data, k):
    kmers = extract_kmers(json_data, k)
    G = nx.MultiDiGraph()
    for kmer in kmers:
        L = kmer[:-1]
        R = kmer[1:]
        G.add_edge(L, R)
    return G


def plot_graph(graph, filename):
    plt.figure(figsize=(8, 8))
    pos = nx.circular_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=300, node_color='lightblue')
    nx.draw_networkx_labels(graph, pos, font_size=8)
    # Draw each edge separately with curvature to visualize multiple edges
    for i, edge in enumerate(graph.edges(keys=True)):
        u, v, k = edge
        nx.draw_networkx_edges(
            graph,
            pos,
            edgelist=[(u, v)],
            connectionstyle=f"arc3,rad={0.15 * (k + 1)}",
            arrows=True
        )
    plt.title("De Bruijn Graph")
    plt.savefig(filename)


def is_connected(graph) -> bool:
    undirected_graph = graph.to_undirected()
    return nx.is_connected(undirected_graph)


def is_valid_graph(graph):
    # Check if the graph is connected
    if not is_connected(graph):
        return False
    unique_degrees = {}
    for node in graph.nodes():
        # Store the degree difference
        diff = graph.out_degree(node) - graph.in_degree(node)
        if diff != 0:
            # Store only nodes with degree difference
            unique_degrees[node] = diff
    # Check if the number of nodes with degree difference is 0 or 2,
    # and if 2, they must be +1 and -1
    if len(unique_degrees) == 0:
        return True
    elif len(unique_degrees) == 2:
        return sorted(unique_degrees.values()) == [-1, 1]
    else:
        return False
