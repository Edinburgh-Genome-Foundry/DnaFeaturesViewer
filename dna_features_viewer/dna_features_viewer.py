import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from .utils import change_luminosity, get_text_box, compute_features_levels


class GraphicFeature:
    """Genetic Feature to be plotted.

    Parameters
    ----------

    start, end
      Coordinates of the feature in the final sequence.

    strand
      Directionality of the feature. can be +1/-1/0 for direct sense,
      anti-sense, or no directionality.

    label
      Short descriptive text associated and plotted with the feature

    color
      Color of the feature, any Matplotlib-compatible format is accepted,
      such as "white", "w", "#ffffff", (1,1,1), etc.

    data
      Any other keyword is kept into the feature.data[] dictionary.
    """

    def __init__(self, start=None, end=None, strand=None,
                 label=None, color="#000080", **data):
        self.start = start
        self.end = end
        self.strand = strand
        self.label = label
        self.color = "#000080" if color is None else color
        self.data = data

    def overlaps_with(self, other):
        """Return True iff the feature's location overlaps with feature `other`
        """
        loc1, loc2 = (self.start, self.end), (other.start, other.end)
        loc1, loc2 = sorted(loc1), sorted(loc2)
        loc1, loc2 = sorted([loc1, loc2], key=lambda loc: loc[0])
        return loc1[1] > loc2[0]

    @property
    def length(self):
        """Return the length of the feature (end-start)"""
        return abs(self.end - self.start)

    @property
    def x_center(self):
        """Return the x-center of the feature, (start+end)/2"""
        return 0.5 * (self.start + self.end)

    @staticmethod
    def from_biopython_feature(feature, color=None, label=None):
        """Create a GraphicalFeature from a Biopython.Feature object."""
        return GraphicFeature(start=feature.location.start,
                              end=feature.location.end,
                              strand=feature.location.strand,
                              color=color, label=label)

    def create_patch(self, level):
        """Create an Arrow Matplotlib patch with the feature's coordinates.

        The Arrow points in the direction of the feature's strand.
        If the feature has no direction (strand==0), the returned patch will
        simply be a rectangle.

        The x-coordinates of the patch are determined by the feature's
        `start` and `end` while the y-coordinates are determined by the `level`
        """
        x1, x2 = self.start, self.end
        if self.strand == -1:
            x1, x2 = x2, x1
        head_length = 5 if self.strand in (-1, 1) else 0
        arrowstyle = mpatches.ArrowStyle.Simple(head_width=14,
                                                tail_width=14,
                                                head_length=head_length)
        return mpatches.FancyArrowPatch([x1, level], [x2, level],
                                        shrinkA=0.0, shrinkB=0.0,
                                        arrowstyle=arrowstyle,
                                        facecolor=self.color)

    def annotate(self, ax, level):
        """Create a Matplotlib Text with the feature's label.

        The x-coordinates of the text are determined by the feature's
        `x_center` while the y-coordinates are determined by the `level`.

        The text is horizontally and vertically centered.

        The Arrow points in the direction of the feature's strand.
        If the feature has no direction (strand==0), the returned patch will
        simply be a rectangle.

        The x-coordinates of the patch are determined by the feature's
        `start` and `end` while the y-coordinates are determined by the `level`
        """
        bg_color = change_luminosity(self.color, luminosity=0.95)
        text = ax.text(
            self.x_center, level, self.label,
            horizontalalignment="center",
            verticalalignment="center",
            bbox=dict(boxstyle="round", fc=bg_color, ec="0.5", lw=1)
        )
        margin = 0.1 * abs(self.end - self.start)
        x1, y1, x2, y2 = get_text_box(text, margin=margin)
        overflowing = (x1 < self.start) or (x2 > self.end)
        return text, overflowing, (x1, x2)

    def plot(self, ax, level=0):
        """Plot the feature's Matplotlib patch on a Matplotlib ax."""
        ax.add_patch(self.create_patch(level=level))

    def __repr__(self):
        return "GF(%(label)s, %(start)d-%(end)d (%(strand)d))" % self.__dict__


class GraphicRecord:
    """Set of Genetic Features of a same DNA sequence, to be plotted together.

    Parameters
    ----------

    sequence_length
      Length of the DNA sequence, in number of nucleotides

    features
      list of GraphicalFeature objects.
    """

    def __init__(self, sequence_length, features):
        self.features = features
        self.sequence_length = sequence_length

    def plot(self, ax=None, fig_width=8):
        """Plot all the features in the same Matplotlib ax

        `fig_width` represents the width in inches of the final figure (if
        no `ax` parameter is attributed).
        """
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
        """Create a new GraphicRecord from a BioPython Record object.

        Parameters
        ----------

        record
          A BioPython Record object

        features_filter
          A function (biofeature->bool) where biofeature is a Biopython Feature
          object found in the record, and bool (True/False) indicates whether
          the feature should be kept (True) of discarded (False).

        fun_color
          A function (biofeature->color) where biofeature is a Biopython
          Feature object found in the record, and color is any
          Matplotlib-compatible color format.

        fun_label
            A function (biofeature->label) where biofeature is a Biopython
            Feature object found in the record, and label is a description
            extracted from the feature.
        """
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
