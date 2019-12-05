"""Implements the method used for deciding which feature goes to which level
when plotting."""

import itertools
import math


class Graph:
    """Minimal implementation of non-directional graphs.

    Parameters
    ----------

    nodes
      A list of objects. They must be hashable
    edges
      A list of the form [(n1,n2), (n3,n4)...] where (n1, n2) represents
      an edge between nodes n1 and n2
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
    - The nodes of the graph are features and they will be colored with a level
    - Two nodes are neighbors if and only if their features's locations overlap
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
    levels = {
        n: None if "fixed_level" not in n.data else n.data["fixed_level"]
        for n in graph.nodes
    }

    def collision(base_level, node):
        """Return whether the node placed at base_level collides with its
        neighbors in the graph."""
        for neighbor in graph.neighbors[node]:
            level = levels[neighbor]
            if level is None:
                continue
            if "nlines" in neighbor.data:
                top = math.ceil(level + 0.5 * neighbor.data["nlines"])
                if level <= base_level < top:
                    return True

                top = math.ceil(base_level + 0.5 * node.data["nlines"])
                if base_level <= level < top:
                    return True
            else:
                if level == base_level:
                    return True
        return False

    for node in sorted(graph.nodes, key=lambda f: -f.length):
        base_level = 0
        while collision(base_level, node):
            base_level += 1
        levels[node] = base_level
    return levels
