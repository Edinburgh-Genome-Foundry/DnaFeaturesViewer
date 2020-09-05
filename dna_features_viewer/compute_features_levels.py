"""Implements the method used for deciding which feature goes to which level
when plotting."""

import itertools
import math


class Graph:
    """Minimal implementation of non-directional graphs.

    Parameters
    ----------

    nodes
      A list of objects. They must be hashable.

    edges
      A list of the form [(n1,n2), (n3,n4)...] where (n1, n2) represents
      an edge between nodes n1 and n2.
    """

    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.neighbors = {n: [] for n in nodes}
        for n1, n2 in edges:
            self.neighbors[n1].append(n2)
            self.neighbors[n2].append(n1)


def compute_features_levels(features):
    """Compute the vertical levels on which the features should be displayed
    in order to avoid collisions.

    `features` must be a list of `dna_features_viewer.GraphicFeature`.

    The method used is basically a graph coloring:
    - The nodes of the graph are features and they will be colored with a level.
    - Two nodes are neighbors if and only if their features's locations overlap.
    - Levels are attributed to nodes iteratively starting with the nodes
      corresponding to the largest features.
    - A node receives the lowest level (starting at 0) that is not already
      the level of one of its neighbors.
    """
    edges = [
        (f1, f2)
        for f1, f2 in itertools.combinations(features, 2)
        if f1.overlaps_with(f2)
    ]
    graph = Graph(features, edges)
    levels = {n: n.data.get("fixed_level", None) for n in graph.nodes}

    def collision(node, level):
        """Return whether the node placed at base_level collides with its
        neighbors in the graph."""
        line_factor = 0.5
        nlines = node.data.get("nlines", 1)
        for neighbor in graph.neighbors[node]:
            neighbor_level = levels[neighbor]
            if neighbor_level is None:
                continue
            neighbor_lines = neighbor.data.get("nlines", 1)
            min_distance = line_factor * (nlines + neighbor_lines)
            if abs(level - neighbor_level) < min_distance:
                return True
        return False

    for node in sorted(graph.nodes, key=lambda f: -f.length):
        if levels[node] is None:
            level = 0
            while collision(node, level):
                level += 0.5
            levels[node] = level
    return levels
