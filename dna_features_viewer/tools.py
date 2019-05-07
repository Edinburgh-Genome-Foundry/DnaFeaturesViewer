"""Useful functions for the library"""

import itertools
import numpy
from matplotlib.colors import colorConverter

def change_luminosity(color, luminosity=None, factor=1):
    """Return a version of the color with different luminosity.

    Parameters
    ----------
    color
      A color in any Matplotlib-compatible format such as "white", "w",
      (1,1,1), "#ffffff", etc.
    luminosity
      A float in 0-1. If provided, the returned color has this level of
      luminosity.
    factor
      Only used if `luminosity` is not set. The luminosity of the new image
      is (1 - (1-L)/(factor + 1)), where L is the current luminosity.
    """
    res = numpy.array(colorConverter.to_rgb(color))
    if luminosity is not None:
        if luminosity == 1:
            return (1, 1, 1)
        l = 1.0*res.sum() / 3.0
        factor = (luminosity - l) / (1.0 - luminosity)
    if factor == -1:
        return numpy.array([1.0, 1.0, 1.0])
    else:
        return (factor + res) / (factor + 1)


def get_text_box(text, margin=0):
    """Return the coordinates of a Matplotlib Text.

    `text` is a Matplotlib text obtained with ax.text().
    This returns `(x1,y1, x2, y2)` where (x1,y1) is the lower left corner
    and (x2, y2) is the upper right corner of the text, in data coordinates.
    If a margin m is supplied, the returned result is (x1-m, y1-m, x2+m, y2+m)
    """
    renderer = text.axes.figure.canvas.get_renderer()
    bbox = text.get_window_extent(renderer)  # bounding box
    bbox_data = bbox
    # bbox_data = bbox.transformed(text.axes.transData.inverted())
    x1, y1, x2, y2 = bbox_data.get_points().flatten()
    return [x1 , y1, x2, y2]

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
    - Two nodes are neighbors iff their features's locations overlap
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
    levels = {n: None for n in graph.nodes}

    def collision(base_level, node):
        """Return True iff the node placed at base_level collides with
        its neighbors in the graph."""
        for neighbor in graph.neighbors[node]:
            level = levels[neighbor]
            if level is None:
                continue
            if 'nlines' in neighbor.data:
                top = numpy.ceil(level + 0.5 * neighbor.data['nlines'])
                if level <= base_level < top:
                    return True
                
                top = numpy.ceil(base_level + 0.5 * node.data['nlines'])
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


def bokeh_feature_patch(self, start, end, strand, width=0.3, level=0, **kw):
    """Return a dict with points coordinates of a Bokeh Feature arrow."""
    hw = width/2.0
    x1, x2 = (start, end) if (strand >= 0) else (end, start)
    if strand >= 0:
        head_base = max(x1, x2 - max(.025*self.sequence_length, .025*(x2-x1)))
    else:
        head_base = min(x1, x2 + max(.025*self.sequence_length, .025*(x1-x2)))
    result = dict(
        xs= [x1, x1, head_base, x2, head_base, x1],
        ys= [e + level for e in [-hw, hw, hw, 0, -hw, -hw]]
    )
    result.update(kw)
    return result