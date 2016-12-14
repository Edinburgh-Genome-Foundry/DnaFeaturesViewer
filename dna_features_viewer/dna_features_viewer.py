import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from .utils import (change_luminosity, get_text_box,
                    compute_features_levels)
import numpy as np
from copy import deepcopy

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
                 label=None, color="#000080", thickness=14, linewidth=1.0,
                 **data):
        self.start = start
        self.end = end
        self.strand = strand
        self.label = label
        self.color = color
        self.data = data
        self.thickness = thickness
        self.linewidth = linewidth

    def split_in_two(self, x_coord=0):
        copy1 = deepcopy(self)
        copy2 = deepcopy(self)
        copy1.end = x_coord
        copy2.start = x_coord + 1
        return copy1, copy2

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
    def from_biopython_feature(feature, **props):
        """Create a GraphicalFeature from a Biopython.Feature object."""
        return GraphicFeature(start=feature.location.start,
                              end=feature.location.end,
                              strand=feature.location.strand,
                              **props)

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

    def plot_feature(self, ax, feature, level, linewidth=1.0):
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
        head_length = 5 if feature.strand in (-1, 1) else 0.001
        arrowstyle = mpatches.ArrowStyle.Simple(head_width=feature.thickness,
                                                tail_width=feature.thickness,
                                                head_length=head_length)
        y = self.feature_level_width * level
        patch = mpatches.FancyArrowPatch([x1, y], [x2, y],
                                         shrinkA=0.0, shrinkB=0.0,
                                         arrowstyle=arrowstyle,
                                         facecolor=feature.color, zorder=0,
                                         linewidth=feature.linewidth)
        ax.add_patch(patch)
        return patch

    def coordinates_in_plot(self, x, level):
        return (x, level * self.feature_level_width)

    def split_overflowing_features_circularly(self):
        """Split the features that overflow over the edge for circular
        constructs (inplace)."""
        new_features = []
        for f in self.features:
            if f.start < 0 < f.end:
                f1, f2 = f.split_in_two(-1)
                f1.start, f1.end = (f1.start + self.sequence_length,
                                    f1.end + self.sequence_length)
                new_features += [f1, f2]
            elif f.start < self.sequence_length < f.end:
                f1, f2 = f.split_in_two(self.sequence_length)
                f2.start, f2.end = (f2.start - self.sequence_length,
                                    f2.end - self.sequence_length)
                new_features += [f1, f2]
            else:
                new_features.append(f)
        self.features = new_features

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

        figure_width = ax.figure.get_size_inches()[0]
        margin = 0.05*figure_width
        x1, y1, x2, y2 = get_text_box(text, margin=margin)
        overflowing = (x1 < feature.start) or (x2 > feature.end)
        return text, overflowing, (x1, x2)

    def plot(self, ax=None, figure_width=8, draw_line=True, with_ruler=True,
             fontsize=11, box_linewidth=1, box_color=None,
             annotate_inline=False, level_offset=0):
        """Plot all the features in the same Matplotlib ax

        `figure_width` represents the width in inches of the final figure (if
        no `ax` parameter is attributed).
        """
        features_levels = compute_features_levels(self.features)
        for f in features_levels:
            features_levels[f] += level_offset
        max_level = (1 if (features_levels == {}) else
                     max(1, max(features_levels.values())))
        auto_figure_height = ax is None
        if ax is None:
            fig, ax = plt.subplots(1, figsize=(figure_width, 2 * max_level))
        self.initialize_ax(ax, draw_line=draw_line, with_ruler=with_ruler)
        overflowing_annotations = []
        for feature, level in features_levels.items():
            self.plot_feature(ax=ax, feature=feature, level=level)
            if feature.label is not None:
                text, overflowing, (x1, x2) = self.annotate_feature(
                    ax=ax, feature=feature, level=level, fontsize=fontsize,
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
                         0 if len(annotations_levels) == 0 else
                         max(annotations_levels.values()),
                         auto_figure_height)
        return ax, annotations_levels

    def finalize_ax(self, ax, features_levels, annotations_levels,
                    auto_figure_height=False):

        ymax = (self.feature_level_width * (features_levels + 1) +
                self.annotation_level_width * (annotations_levels + 1))
        ax.set_ylim(-1, ymax)
        if auto_figure_height:
            figure_width = ax.figure.get_size_inches()[0]
            ax.figure.set_size_inches(figure_width, 1.5 + 0.37 * ymax)

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
            return mpatches.Wedge._recompute_path(self)

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
                 feature_level_width=0.2, annotation_level_width=0.25):

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

        ax.set_xlim(-self.radius, self.radius)
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
        ratio = 1.0 * (ymax - ymin) / (xmax - xmin)

        if auto_figure_height:
            figure_width = ax.figure.get_size_inches()[0]
            ax.figure.set_size_inches(figure_width, figure_width * ratio)

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

class BiopythonTranslator:

    def __init__(self, features_filters=(), features_properties=None):
        self.features_filters = features_filters
        self.features_properties = features_properties

    @staticmethod
    def compute_feature_color(feature):
        return feature.qualifiers.get("color", "#7245dc")

    def compute_filtered_features(self, features):
        return [f for f in features
                if all([fl(f) for fl in self.features_filters])]

    @staticmethod
    def compute_feature_label(feature):
        """Gets the 'label' of the feature."""
        result = feature.type
        for key in ["label", "source", "locus_tag", "note"]:
            if key in feature.qualifiers:
                result = feature.qualifiers[key]
        if isinstance(result, list):
            return "|".join(result)
        else:
            return result

    def translate_feature(self, feature):
        properties = dict(label=self.compute_feature_label(feature),
                          color=self.compute_feature_color(feature))
        if self.features_properties is not None:
            other_properties = self.features_properties(feature)
        else:
            other_properties = {}
        properties.update(other_properties)
        return GraphicFeature(start=feature.location.start,
                              end=feature.location.end,
                              strand=feature.location.strand,
                              **properties)

    def translate_record(self, record, grecord_class=None):
        """Create a new GraphicRecord from a BioPython Record object.

        Parameters
        ----------

        record
          A BioPython Record object
        """
        if grecord_class is None:
            grecord_class = GraphicRecord
        return grecord_class(sequence_length=len(record.seq), features=[
            self.translate_feature(feature)
            for feature in self.compute_filtered_features(record.features)
            if feature.location is not None
        ])
