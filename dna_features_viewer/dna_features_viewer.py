import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import itertools as itt
from matplotlib.colors import colorConverter
import numpy as np

def lighten(color, factor=1, luminosity=None):
    res = np.array(colorConverter.to_rgb(color))
    if luminosity is not None:
        if luminosity == 1:
            return (1, 1, 1)
        l = res.sum() / 3.0
        factor = (luminosity - l) / (1.0 - luminosity)
    return (factor + res) / (factor + 1)


class Graph:

    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.neighbors = {n: [] for n in nodes}
        for n1, n2 in edges:
            self.neighbors[n1].append(n2)
            self.neighbors[n2].append(n1)


def compute_features_levels(features):
    edges = [
        (f1, f2)
        for f1, f2 in itt.combinations(features, 2)
        if f1.overlaps_with(f2)
    ]
    graph = Graph(features, edges)
    levels = {
        n: None
        for n in graph.nodes
    }
    for node in sorted(graph.nodes, key=lambda f: -f.length):
        level = 0
        while any([levels[n] == level
                   for n in graph.neighbors[node]]):
            level += 1
        levels[node] = level
    return levels


def get_text_box(text, margin=0):
    renderer = text.axes.figure.canvas.get_renderer()
    bbox = text.get_window_extent(renderer)
    bbox_data = bbox.transformed(text.axes.transData.inverted())
    x1, y1, x2, y2 = bbox_data.get_points().flatten()
    return [x1 - margin, y1 - margin, x2 + margin, y2 + margin]


class GraphicFeature:

    def __init__(self, start=None, end=None, strand=None,
                 label=None, color="#000080", **data):
        self.start = start
        self.end = end
        self.strand = strand
        self.label = label
        self.color = "#000080" if color is None else color
        self.data = data

    def overlaps_with(self, other):
        loc1, loc2 = (self.start, self.end), (other.start, other.end)
        loc1, loc2 = sorted(loc1), sorted(loc2)
        loc1, loc2 = sorted([loc1, loc2], key=lambda loc: loc[0])
        return loc1[1] > loc2[0]

    @property
    def length(self):
        return abs(self.end - self.start)

    @property
    def x_center(self):
        return 0.5 * (self.start + self.end)

    @staticmethod
    def from_biopython_feature(feature, color=None, label=None):
        return GraphicFeature(start=feature.location.start,
                              end=feature.location.end,
                              strand=feature.location.strand,
                              color=color, label=label)

    def create_patch(self, level):
        x1, x2 = self.start, self.end
        if self.strand == -1:
            x1, x2 = x2, x1
        head_length = 5 if self.strand in (-1, 1) else 0
        return mpatches.FancyArrowPatch(
            [x1, level], [x2, level],
            arrowstyle=mpatches.ArrowStyle.Simple(
                head_width=14, tail_width=14, head_length=head_length
            ), shrinkA=0.0, shrinkB=0.0,
            facecolor=self.color
        )

    def annotate(self, ax, level):
        text = ax.text(self.x_center, level, self.label,
                       horizontalalignment="center",
                       verticalalignment="center",
                       bbox=dict(boxstyle="round",
                                 fc=lighten(self.color, luminosity=0.95),
                                 ec="0.5", lw=1)
                       )
        x1, y1, x2, y2 = get_text_box(text)
        padding = 0.1 * abs(self.end - self.start)
        x1, x2 = x1 - padding, x2 + padding
        overflowing = (x1 < self.start) or (x2 > self.end)
        return text, overflowing, (x1, x2)

    def plot(self, ax, level=0):
        ax.add_patch(self.create_patch(level=level))

    def __repr__(self):
        return "GF(%(label)s, %(start)d-%(end)d (%(strand)d))" % self.__dict__


class GraphicRecord:

    def __init__(self, sequence_length, features):
        self.features = features
        self.sequence_length = sequence_length

    def plot(self, ax=None, fig_width=8):
        levels = compute_features_levels(self.features)
        max_level = max(levels.values())
        auto_figure_height = ax is None
        if ax is None:
            fig, ax = plt.subplots(1, figsize=(fig_width, 2 * max_level))
            ax.set_ylim(-1, max_level + 1),
            ax.set_xlim(0, self.sequence_length)
        ax.axis("off")
        ax.set_xlim(0, self.sequence_length)
        ax.set_ylim(-1, max_level + 1)
        overflowing_annotations = []
        for feature, level in levels.items():
            feature.plot(ax=ax, level=level)
            text, overflowing, (x1, x2) = feature.annotate(ax=ax, level=level)
            if overflowing:
                ft = GraphicFeature(start=x1, end=x2, feature=feature,
                                    text=text, feature_level=level)
                overflowing_annotations.append(ft)
        levels = compute_features_levels(overflowing_annotations)
        max_y = max_level + 1
        for feature, level in levels.items():
            text = feature.data["text"]
            ft = feature.data["feature"]
            x, y = text.get_position()
            new_y = max_level + 1 + level
            max_y = max(max_y, new_y)
            text.set_position((x, new_y))
            xx = [x, feature.data["feature"].x_center]
            yy = [new_y, feature.data["feature_level"]]
            ax.plot(xx, yy, c="#cccccc", lw=1, zorder=-1000,)
        ax.set_ylim(-1, max_y + 1)
        if auto_figure_height:
            ax.figure.set_size_inches(fig_width, 1.5 + 0.37 * max_y)
        return ax, max_y

    @staticmethod
    def from_biopython_record(record, features_filter=None,
                              fun_color=None, fun_label=None):
        if fun_color is None:
            fun_color = lambda a: None
        if fun_label is None:
            fun_label = lambda f: f.type
        if features_filter is None:
            features_filter = lambda a: True

        features = [
            GraphicFeature.from_biopython_feature(
                feature,
                color=fun_color(feature),
                label=fun_label(feature)
            )
            for feature in record.features
            if features_filter(feature)
        ]
        return GraphicRecord(sequence_length=len(record.seq),
                             features=features)
