import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from .utils import change_luminosity, get_text_box, compute_features_levels

try:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.Alphabet import DNAAlphabet
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


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
    feature_type = "feature"

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

    def __repr__(self):
        return (("GF(%(label)s, %(start)d-%(end)d " % self.__dict__) +
                (")" if self.strand is None else "(%d))" % self.strand))


class GraphicRecord:
    """Set of Genetic Features of a same DNA sequence, to be plotted together.

    Parameters
    ----------

    sequence_length
      Length of the DNA sequence, in number of nucleotides

    features
      list of GraphicalFeature objects.
    """

    def __init__(self, sequence_length, features, feature_level_width=1,
                 annotation_level_width=1):
        self.features = features
        self.sequence_length = sequence_length
        self.feature_level_width = feature_level_width
        self.annotation_level_width = annotation_level_width

    def initialize_ax(self, ax, draw_line, with_ruler):

        if draw_line:
            ax.axhline(0, zorder=-1000, c="k")

        if with_ruler:  # only display the xaxis ticks
            ax.set_frame_on(False)
            ax.yaxis.set_visible(False)
            ax.xaxis.tick_bottom()
        else:  # don't display anything
            ax.axis("off")

        ax.set_xlim(0, self.sequence_length)

    def plot_feature(self, ax, feature, level):
        """Create an Arrow Matplotlib patch with the feature's coordinates.

        The Arrow points in the direction of the feature's strand.
        If the feature has no direction (strand==0), the returned patch will
        simply be a rectangle.

        The x-coordinates of the patch are determined by the feature's
        `start` and `end` while the y-coordinates are determined by the `level`
        """
        x1, x2 = feature.start, feature.end
        if feature.strand == -1:
            x1, x2 = x2, x1
        head_length = 5 if feature.strand in (-1, 1) else 0
        arrowstyle = mpatches.ArrowStyle.Simple(head_width=14,
                                                tail_width=14,
                                                head_length=head_length)
        patch = mpatches.FancyArrowPatch([x1, level], [x2, level],
                                         shrinkA=0.0, shrinkB=0.0,
                                         arrowstyle=arrowstyle,
                                         facecolor=feature.color, zorder=0)
        ax.add_patch(patch)
        return patch

    def coordinates_in_plot(self, x, level):
        return (x, level * self.feature_level_width)

    def annotate_feature(self, ax, feature, level, fontsize=11,
                         box_linewidth=1, box_color=None):
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
        bg_color = change_luminosity(feature.color, luminosity=0.95)
        x, y = self.coordinates_in_plot(feature.x_center, level)
        text = ax.text(
            x, y, feature.label,
            horizontalalignment="center",
            verticalalignment="center",
            bbox=dict(boxstyle="round", fc=bg_color if box_color is None
                      else box_color, ec="0.5",
                      lw=box_linewidth),
            fontsize=fontsize,
            zorder=2
        )
        margin = 0.1 * abs(feature.end - feature.start)
        x1, y1, x2, y2 = get_text_box(text, margin=margin)
        overflowing = (x1 < feature.start) or (x2 > feature.end)
        return text, overflowing, (x1, x2)

    def plot(self, ax=None, fig_width=8, draw_line=True, with_ruler=True,
             fontsize=11, box_linewidth=1, box_color=None,
             annotate_inline=False):
        """Plot all the features in the same Matplotlib ax

        `fig_width` represents the width in inches of the final figure (if
        no `ax` parameter is attributed).
        """
        features_levels = compute_features_levels(self.features)
        max_level = (1 if (features_levels == {}) else
                     max(1, max(features_levels.values())))
        auto_figure_height = ax is None
        if ax is None:
            fig, ax = plt.subplots(1, figsize=(fig_width, 2 * max_level))
        self.initialize_ax(ax, draw_line=draw_line, with_ruler=with_ruler)
        overflowing_annotations = []
        for feature, level in features_levels.items():
            self.plot_feature(ax=ax, feature=feature, level=level)
            if feature.label is not None:
                text, overflowing, (x1, x2) = self.annotate_feature(
                    ax=ax, feature=feature, level=level,
                    box_linewidth=box_linewidth, box_color=box_color
                )
                if overflowing or not annotate_inline:
                    overflowing_annotations.append(GraphicFeature(
                        start=x1, end=x2, feature=feature,
                        text=text, feature_level=level
                    ))

        annotations_levels = compute_features_levels(overflowing_annotations)
        for feature, level in annotations_levels.items():
            text = feature.data["text"]
            x, y = text.get_position()
            new_y = ((max_level + 1) * self.feature_level_width +
                     (level) * self.annotation_level_width)
            text.set_position((x, new_y))
            fx, fy = self.coordinates_in_plot(feature.data["feature"].x_center,
                                              feature.data["feature_level"])
            ax.plot([x, fx], [new_y, fy], c="k", lw=0.5, zorder=1)

        self.finalize_ax(ax, max(features_levels.values()),
                         max(annotations_levels.values()),
                         auto_figure_height)
        return ax

    def finalize_ax(self, ax, features_levels, annotations_levels,
                    auto_figure_height=False):

        ymax = (self.feature_level_width * (features_levels + 1) +
                self.annotation_level_width * (annotations_levels + 1))
        ax.set_ylim(-1, ymax)
        if auto_figure_height:
            fig_width = ax.figure.get_size_inches()[0]
            ax.figure.set_size_inches(fig_width, 1.5 + 0.37 * ymax)

    @classmethod
    def from_biopython_record(cls, record, features_filter=None,
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
            if feature.location is not None
            if features_filter(feature)
        ]
        return cls(sequence_length=len(record.seq),
                             features=features)

    def to_biopython_record(self):
        """
        Example
        -------
        from Bio import SeqIO
        gr_record = GraphicRecord(features=features, sequence_length=len(seq),
                                  sequence=seq)
        bio_record = gr_record.to_biopython_record()
        with open("example.gb", "w+") as f:
            SeqIO.write(record, f, "genbank")
        """
        if not BIOPYTHON_AVAILABLE:
            raise ImportError(".to_biopython_record requires Biopython")
        features = [
            SeqFeature(FeatureLocation(f.start, f.end, f.strand),
                       type=f.feature_type, qualifiers={"label": f.label})
            for f in self.features
        ]
        sequence = Seq(self.data["sequence"], alphabet=DNAAlphabet())
        return SeqRecord(sequence=sequence, features=features)


import matplotlib.patches as mpatches
from dna_features_viewer import GraphicRecord
import matplotlib.pyplot as plt
import numpy as np


class ArrowWedge(mpatches.Wedge):

    def __init__(self, center, radius, theta1, theta2, width, direction=+1,
                 **kwargs):

        self.direction = direction
        self.radius = radius
        mpatches.Wedge.__init__(self, center, radius,
                                theta1, theta2, width, **kwargs)
        self._recompute_path()

    def _recompute_path(self):

        if not self.direction in [-1, +1]:
            return Wedge._recompute_path(self)

        theta1, theta2 = self.theta1, self.theta2
        arrow_angle = min(5, abs(theta2 - theta1) / 2)
        normalized_arrow_width = self.width / 2.0 / self.radius
        if self.direction == +1:
            angle_start_arrow = theta1 + arrow_angle
            arc = mpatches.Path.arc(angle_start_arrow, theta2)
            outer_arc = arc.vertices[::-1] * (1 + normalized_arrow_width)
            inner_arc = arc.vertices * (1 - normalized_arrow_width)
            arrow_vertices = [
                outer_arc[-1],
                np.array([np.cos(np.deg2rad(theta1)),
                          np.sin(np.deg2rad(theta1))]),
                inner_arc[0]
            ]
        else:
            angle_start_arrow = theta2 - arrow_angle
            arc = mpatches.Path.arc(theta1, angle_start_arrow)
            outer_arc = arc.vertices * \
                (self.radius + self.width / 2.0) / self.radius
            inner_arc = arc.vertices[
                ::-1] * (self.radius - self.width / 2.0) / self.radius
            arrow_vertices = [
                outer_arc[-1],
                np.array([np.cos(np.deg2rad(theta2)),
                          np.sin(np.deg2rad(theta2))]),
                inner_arc[0]
            ]
        p = np.vstack([outer_arc, arrow_vertices, inner_arc])

        path_vertices = np.vstack([p, inner_arc[-1, :], (0, 0)])

        path_codes = np.hstack([arc.codes,
                                4 * [mpatches.Path.LINETO],
                                arc.codes[1:],
                                mpatches.Path.LINETO,
                                mpatches.Path.CLOSEPOLY])
        path_codes[len(arc.codes)] = mpatches.Path.LINETO

        # Shift and scale the wedge to the final location.
        path_vertices *= self.r
        path_vertices += np.asarray(self.center)
        self._path = mpatches.Path(path_vertices, path_codes)


class CircularGraphicRecord(GraphicRecord):

    def __init__(self, sequence_length, features, top_position=0,
                 feature_level_width=0.2, annotation_level_width=0.2):

        self.radius = 1.0
        self.sequence_length = sequence_length
        self.features = features
        self.top_position = top_position
        self.angle_rotation = 2 * np.pi * top_position / sequence_length
        self.feature_level_width = feature_level_width
        self.annotation_level_width = annotation_level_width

    def initialize_ax(self, ax, draw_line, with_ruler):

        if draw_line:
            circle = mpatches.Circle(
                (0, -self.radius), self.radius, facecolor='none')
            ax.add_patch(circle)
        ax.axis("off")
        if with_ruler:  # only display the xaxis ticks
            ax.set_frame_on(False)
            ax.yaxis.set_visible(False)
            ax.xaxis.tick_bottom()
        else:  # don't display anything
            ax.axis("off")

        ax.set_xlim(-self.radius,self.radius)
        ax.set_aspect("equal")

    def finalize_ax(self, ax, features_levels, annotation_levels,
                    auto_figure_height=False):
        ymin = (-2 * self.radius -
                self.feature_level_width * (features_levels + 1))
        ymax = (self.feature_level_width * (features_levels + 1) +
                (annotation_levels + 1) * self.annotation_level_width)
        xmin = -self.radius - self.feature_level_width * (features_levels + 1)
        xmax = -xmin
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ratio = 1.0*(ymax - ymin)/(xmax - xmin)

        if auto_figure_height:
            fig_width = ax.figure.get_size_inches()[0]
            ax.figure.set_size_inches(fig_width, fig_width *ratio)

    def plot_feature(self, ax, feature, level):
        a_start = self.position_to_angle(feature.start)
        a_end = self.position_to_angle(feature.end)
        a_start, a_end = sorted([a_start, a_end])
        r = self.radius + level * self.feature_level_width
        patch = ArrowWedge((0, -self.radius), r, a_start, a_end,
                           0.7 * self.feature_level_width,
                           direction=feature.strand,
                           facecolor=feature.color, zorder=1)
        ax.add_patch(patch)

    def position_to_angle(self, position):
        a = 360.0 * (position - self.top_position) / self.sequence_length
        return 90 - a

    def coordinates_in_plot(self, x, level):
        r = self.radius + level * self.feature_level_width
        angle = self.position_to_angle(x)
        rad_angle = np.deg2rad(angle)
        return np.array([r * np.cos(rad_angle),
                         r * np.sin(rad_angle) - self.radius])
