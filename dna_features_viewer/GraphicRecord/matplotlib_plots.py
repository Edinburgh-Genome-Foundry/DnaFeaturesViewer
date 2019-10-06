"""Useful functions for the library"""

import colorsys

import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

from ..biotools import extract_graphical_translation
from ..compute_features_levels import compute_features_levels
from ..GraphicFeature import GraphicFeature
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
        mean = res.sum() / 3.0
        factor = (luminosity - mean) / (1.0 - luminosity)
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
    x1, y1, x2, y2 = bbox_data.get_points().flatten()
    bbox_data = bbox.transformed(text.axes.transData.inverted())
    x1, _, x2, _ = bbox_data.get_points().flatten()
    return [x1, y1, x2, y2]


class MatplotlibPlottableMixin:
    """Class mixin for matplotlib-related methods."""

    def initialize_ax(self, ax, draw_line, with_ruler, ruler_color=None):
        """Initialize the ax: remove axis, draw a horizontal line, etc.

        Parameters
        ----------

        draw_line
          True/False to draw the horizontal line or not.

        with_ruler
          True/False to draw the indices indicators along the line.

        """
        ruler_color = ruler_color or self.default_ruler_color
        start, end = self.span
        plot_start, plot_end = start - 0.8, end - 0.2
        if draw_line:
            ax.plot([plot_start, plot_end], [0, 0], zorder=-1000, c="k")

        if with_ruler:  # only display the xaxis ticks
            ax.set_frame_on(False)
            ax.yaxis.set_visible(False)
            ax.xaxis.tick_bottom()
            if ruler_color is not None:
                ax.tick_params(axis="x", colors=ruler_color)
        else:  # don't display anything
            ax.axis("off")

        ax.set_xlim(plot_start, plot_end)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    def finalize_ax(
        self,
        ax,
        features_levels,
        annotations_max_level,
        auto_figure_height=False,
        ideal_yspan=None,
    ):
        """Redefine y-bounds and figure height.
        """

        # Compute the "natural" ymax
        annotation_height = self.determine_annotation_height(None)
        features_ymax = self.feature_level_height * (features_levels + 1)
        annotations_ymax = annotation_height * annotations_max_level
        ymax = features_ymax + annotations_ymax
        # ymax = self.feature_level_height * (
        #     features_levels + 1
        # ) + annotation_height * (
        #     annotations_max_level + (annotations_max_level > 0)
        # )
        ymin = -1

        # ymax could be even bigger if a "ideal_yspan" has been set.
        if (ideal_yspan is not None) and not (auto_figure_height):
            ymax = max(ideal_yspan + ymin, ymax)
        ax.set_ylim(ymin, ymax)
        if auto_figure_height:
            figure_width = ax.figure.get_size_inches()[0]
            ax.figure.set_size_inches(figure_width, 1 + 0.4 * ymax)
        if self.plots_indexing == "genbank":
            ax.set_xticklabels([int(i + 1) for i in ax.get_xticks()])
        return ideal_yspan / (ymax - ymin)

    def plot_feature(self, ax, feature, level, linewidth=1.0):
        """Create an Arrow Matplotlib patch with the feature's coordinates.

        The Arrow points in the direction of the feature's strand.
        If the feature has no direction (strand==0), the returned patch will
        simply be a rectangle.

        The x-coordinates of the patch are determined by the feature's
        `start` and `end` while the y-coordinates are determined by the `level`
        """
        x1, x2 = feature.start, feature.end
        if feature.open_left:
            x1 -= 1
        if feature.open_right:
            x2 += 1
        if feature.strand == -1:
            x1, x2 = x2, x1
        x1, x2 = x1 - 0.5, x2 - 0.5

        is_undirected = feature.strand not in (-1, 1)
        head_is_cut = (feature.strand == 1 and feature.open_right) or (
            feature.strand == -1 and feature.open_left
        )
        head_length = (
            0.001
            if (is_undirected or head_is_cut)
            else max(0.6 * feature.thickness, 5)
        )

        arrowstyle = mpatches.ArrowStyle.Simple(
            head_width=feature.thickness,
            tail_width=feature.thickness,
            head_length=head_length,
        )
        y = self.feature_level_height * level
        patch = mpatches.FancyArrowPatch(
            [x1, y],
            [x2, y],
            shrinkA=0.0,
            shrinkB=0.0,
            arrowstyle=arrowstyle,
            facecolor=feature.color,
            zorder=0,
            edgecolor=feature.linecolor,
            linewidth=feature.linewidth,
        )
        ax.add_patch(patch)
        return patch

    def autoselect_label_color(self, background_color):
        """Autselect a color for the label font.

        In the current method the label will be black on clear backrgounds,
        and white on dark backgrounds.
        """
        r, g, b = colorConverter.to_rgb(background_color)
        luminosity = 0.299 * r + 0.587 * g + 0.114 * b
        return "black" if (luminosity > 0.4) else "white"

    def annotate_feature(
        self,
        ax,
        feature,
        level,
        inline=False,
        max_label_length=50,
        max_line_length=30,
        padding=0,
    ):
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
        label = feature.label
        if not inline:
            label = self._format_label(
                label,
                max_label_length=max_label_length,
                max_line_length=max_line_length,
            )
        nlines = len(label.split("\n"))
        fontdict = dict(**feature.fontdict)
        if "family" not in fontdict and (self.default_font_family is not None):
            fontdict["family"] = self.default_font_family
        if inline and ("color" not in fontdict):
            color = self.autoselect_label_color(background_color=feature.color)
            fontdict["color"] = color
        box_color = feature.box_color
        if box_color == "auto":
            box_color = self.default_box_color
        bbox = None
        if (box_color is not None) and not inline:
            bbox = dict(
                boxstyle="round",
                fc=bg_color if box_color == "auto" else box_color,
                ec="0.5",
                lw=feature.box_linewidth,
            )
        text = ax.text(
            x,
            y,
            label,
            horizontalalignment="center",
            verticalalignment="center",
            bbox=bbox,
            fontdict=fontdict,
            zorder=2,
        )
        x1, y1, x2, y2 = get_text_box(text)
        x1 -= padding
        x2 += padding
        overflowing = (x1 < feature.start) or (x2 > feature.end)
        return text, overflowing, nlines, (x1, x2), (y2 - y1)

    def position_annotation(self, feature, ax, level, annotate_inline):
        padding = self.compute_padding(ax)
        if annotate_inline:
            text, overflowing, lines, (x1, x2), height = self.annotate_feature(
                ax=ax,
                feature=feature,
                level=level,
                inline=True,
                padding=padding,
            )
            
            if overflowing:
                text.remove()
                text, _, lines, (x1, x2), height = self.annotate_feature(
                    ax=ax,
                    feature=feature,
                    level=level,
                    inline=False,
                    padding=padding,
                )
            return text, overflowing, lines, (x1, x2), height
        else:
            return self.annotate_feature(
                ax=ax, feature=feature, level=level, padding=padding
            )

    def plot(
        self,
        ax=None,
        figure_width=8,
        draw_line=True,
        with_ruler=True,
        ruler_color=None,
        plot_sequence=False,
        annotate_inline=True,
        max_label_length=50,
        max_line_length=30,
        level_offset=0,
        x_lim=None,
        figure_height=None,
    ):
        """Plot all the features in the same Matplotlib ax

        Parameters
        ----------

        ax
          The Matplotlib ax on which to plot the graphic record. If None is
          provided, a new figure and ax is generated, the ax is returned at
          the end.

        figure_width
          Width of the figure (only if no ax was provided and a new figure is
          created) in inches.

        draw_line
          If True, a base line representing the sequence will be drawn.

        with_ruler
          If true, the sequence indices will be indicated at regular intervals.

        plot_sequence
          If True and the graphic record has a "sequence" attribute set, the
          sequence will be displayed below the base line.

        annotate_inline
          If true, some feature labels will be displayed inside their
          corresponding feature if there is sufficient space.

        level_offset
          All features and annotations will be pushed up by "level_offset". Can
          be useful when plotting several sets of features successively on a
          same ax.

        x_lim
          Horizontal axis limits to be set at the end.
        """
        features_levels = compute_features_levels(self.features)
        for f in features_levels:
            features_levels[f] += level_offset
        max_level = (
            1
            if (features_levels == {})
            else max(1, max(features_levels.values()))
        )
        auto_figure_height = (ax is None) and (figure_height is None)
        if ax is None:
            height = figure_height or max_level
            fig, ax = plt.subplots(1, figsize=(figure_width, height))
        self.initialize_ax(ax, draw_line=draw_line, with_ruler=with_ruler)
        if x_lim is not None:
            ax.set_xlim(*x_lim)
        overflowing_annotations = []
        renderer = ax.figure.canvas.get_renderer()
        bbox = ax.get_window_extent(renderer)
        ax_height = bbox.height
        ideal_yspan = 0
        for feature, level in features_levels.items():
            self.plot_feature(ax=ax, feature=feature, level=level)
            if feature.label is None:
                continue
            text, overflowing, nlines, (
                x1,
                x2,
            ), height = self.position_annotation(
                feature, ax, level, annotate_inline
            )
            # print (height)
            # nlines = len(feature.label.split("\n"))
            line_height = height / nlines
            # Hey look, a magic number!
            feature_ideal_span = 0.4 * (ax_height / line_height)
            ideal_yspan = max(ideal_yspan, feature_ideal_span)
            # print ("ideal_yspan", ideal_yspan, line_height, ax_height)
            if overflowing or not annotate_inline:
                overflowing_annotations.append(
                    GraphicFeature(
                        start=x1,
                        end=x2,
                        feature=feature,
                        text=text,
                        feature_level=level,
                        nlines=nlines,
                    )
                )
        annotations_levels = compute_features_levels(overflowing_annotations)
        max_annotations_level = max([0] + list(annotations_levels.values()))
        annotation_height = self.determine_annotation_height(max_level)
        labels_data = {}
        for feature, level in annotations_levels.items():
            text = feature.data["text"]
            x, y = text.get_position()
            new_y = (max_level + 1) * self.feature_level_height + (
                level
            ) * annotation_height
            text.set_position((x, new_y))
            fx, fy = self.coordinates_in_plot(
                feature.data["feature"].x_center, feature.data["feature_level"]
            )
            ax.plot([x, fx], [new_y, fy], c="k", lw=0.5, zorder=1)
            labels_data[feature.data["feature"]] = dict(
                feature_y=fy, annotation_y=new_y
            )

        if plot_sequence:
            self.plot_sequence(ax)

        self.finalize_ax(
            ax=ax,
            features_levels=max([1] + list(features_levels.values())),
            annotations_max_level=max_annotations_level,
            auto_figure_height=auto_figure_height,
            ideal_yspan=ideal_yspan,
        )
        return ax, (features_levels, labels_data)

    def plot_sequence(
        self,
        ax,
        location=None,
        y_offset=1,
        fontdict=None,
        background=("#f7fbff", "#fffcf0"),
    ):
        """Plot a sequence of nucleotides at the bottom of the plot.

        Parameters
        ----------

        ax
          Which axes the translation should be plotted to

        location
          location of the segment to translate, either (start, end) or
          (start, end, strand)

        y_offset
          Number of text levels under the plot's base line where to draw the
          nucleotides. Should be 1 if the nucleotide sequence is to be plotted
          directly under the main line.

        fontdict
          Matplotlib fontdict for the text, e.g.
          ``{'size': 11, 'weight':'bold'}``

        background
          tuple (color1, color2) of alternate colors to plot behind each
          nucleotide position to guide vision. Leave to None for no background.

        translation
          Sequence of amino acids either as a string ``'MAKG...'`` or as a list
          ``['Met', 'Ala', ...]``
        """
        if self.sequence is None:
            raise ValueError("No sequence in the graphic record")
        if location is None:
            location = self.span
        location_start, location_end = location
        fontdict = dict(size=11)
        fontdict.update(fontdict or {})
        for i, nucleotide in enumerate(self.sequence):
            index = i + location_start
            if location_start <= index <= location_end:
                ax.text(
                    index,
                    -0.7 * self.feature_level_height * y_offset,
                    nucleotide,
                    ha="center",
                    va="center",
                    fontdict=fontdict,
                )
        if background is not None:
            for i in range(location_start, location_end):
                ax.fill_between(
                    [i - 0.5, i + 0.5],
                    y1=1000,
                    y2=-1000,
                    zorder=-2000,
                    facecolor=background[i % 2],
                )
        ymin = ax.get_ylim()[0]
        ax.set_ylim(bottom=min(ymin, -y_offset * self.feature_level_height))

    def plot_translation(
        self,
        ax,
        location=None,
        y_offset=2,
        fontdict=None,
        background=("#f5fff0", "#fff7fd"),
        translation=None,
        long_form_translation=True,
    ):
        """Plot a sequence of amino-acids at the bottom of the plot.

        Parameters
        ----------

        ax
          Which axes the translation should be plotted to

        location
          location of the segment to translate (start, end)

        y_offset
          Number of text levels under the plot's base line where to draw the
          amino acid names. Should be 2 if the nucleotide sequence is also
          plotted at level 1.

        fontdict
          Matplotlib fontdict for the text, e.g.
          ``{'size': 11, 'weight':'bold'}``

        background
          tuple (color1, color2) of alternate colors to plot behind each
          amino acid position to guide vision. Leave to None for no background.

        translation
          Sequence of amino acids either as a string ``'MAKG...'`` or as a list
          ``['Met', 'Ala', ...]``
        """
        start, end = location[0], location[1]
        strand = location[2] if (len(location) == 3) else 1
        s, e = self.span
        start = max(start, s + ((start - s) % 3))
        end = min(end, e - ((end - e) % 3))
        if translation is None:
            new_loc = start - self.first_index, end - self.first_index, strand
            translation = extract_graphical_translation(
                self.sequence,
                location=new_loc,
                long_form=long_form_translation,
            )
        texts = [
            ((start + 3 * i, start + 3 * (i + 1)), aa)
            for i, aa in enumerate(translation)
        ]

        y = -0.7 * y_offset * self.feature_level_height
        ymin = ax.get_ylim()[0]
        ax.set_ylim(bottom=min(ymin, -y_offset * self.feature_level_height))
        fontdict = fontdict or {}
        for i, ((start, end), text) in enumerate(texts):
            ax.text(
                0.5 * (start + end - 1),
                y,
                text,
                ha="center",
                va="center",
                fontdict=fontdict,
            )
            if background:
                ax.fill_between(
                    [start - 0.5, end - 0.5],
                    y1=1000,
                    y2=-1000,
                    zorder=-1000,
                    facecolor=background[i % 2],
                )
