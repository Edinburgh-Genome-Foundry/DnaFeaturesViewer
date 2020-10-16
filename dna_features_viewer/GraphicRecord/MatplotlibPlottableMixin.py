"""Useful functions for the library"""

import colorsys

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
import matplotlib.ticker as ticker

from ..compute_features_levels import compute_features_levels
from ..GraphicFeature import GraphicFeature
from matplotlib.colors import colorConverter
from .MultilinePlottableMixin import MultilinePlottableMixin
from .SequenceAndTranslationMixin import SequenceAndTranslationMixin


class MatplotlibPlottableMixin(MultilinePlottableMixin, SequenceAndTranslationMixin):
    """Class mixin for matplotlib-related methods."""

    default_elevate_outline_annotations = False
    default_strand_in_label_threshold = None

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
        if self.first_index != 0:
            ax.ticklabel_format(useOffset=False, style="plain")
        fmt = lambda x, p: "{:,}".format(int(x))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(fmt))
        if self.ticks_resolution == "auto":
            ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        else:
            locator = ticker.MultipleLocator(self.ticks_resolution)
            ax.xaxis.set_major_locator(locator)

    def finalize_ax(
        self,
        ax,
        features_levels,
        annotations_max_level,
        auto_figure_height=False,
        ideal_yspan=None,
        annotations_are_elevated=True,
    ):
        """Prettify the figure with some last changes.

        Changes include redefining y-bounds and figure height.

        Parameters
        ----------

        ax
          ax on which the record was plotted.

        features_levels

        annotations_max_level
          Number indicating to the method the maximum height for an
          annotation, so the method can set ymax accordingly.

        auto_figure_height
          If true, the figure height will be automatically re-set to a nice
          value (counting ~0.4 inch per level in the figure).

        ideal_yspan
          if provided, can help the method select a better ymax to make sure
          all constraints fit.
        """

        # Compute the "natural" ymax
        annotation_height = self.determine_annotation_height(None)
        features_ymax = self.feature_level_height * (features_levels + 1)
        annotations_ymax = annotation_height * annotations_max_level
        if annotations_are_elevated:
            ymax = features_ymax + annotations_ymax
        else:
            ymax = max(features_ymax, annotations_ymax) + 1
        ymin = min(ax.get_ylim()[0], -0.5)

        # ymax could be even bigger if a "ideal_yspan" has been set.
        if (ideal_yspan is not None) and not (auto_figure_height):
            ymax = max(ideal_yspan + ymin, ymax)
        ax.set_ylim(ymin, ymax)
        if auto_figure_height:
            figure_width = ax.figure.get_size_inches()[0]
            ax.figure.set_size_inches(figure_width, 1 + 0.4 * ymax)
        ax.set_xticks(
            [
                t
                for t in ax.get_xticks()
                if t <= self.last_index and t >= self.first_index
            ]
        )
        if self.plots_indexing == "genbank":
            ax.set_xticklabels([int(i + 1) for i in ax.get_xticks()])
        return ideal_yspan / (ymax - ymin)

    @staticmethod
    def _get_ax_width(ax, unit="inch"):
        """Return the ax's width in 'inches' or 'pixel'."""
        transform = ax.figure.dpi_scale_trans.inverted()
        bbox = ax.get_window_extent().transformed(transform)
        width = bbox.width
        if unit == "pixel":
            width *= ax.figure.dpi
        return width

    def plot_feature(self, ax, feature, level, linewidth=1.0):
        """Create an Arrow Matplotlib patch with the feature's coordinates.

        The Arrow points in the direction of the feature's strand.
        If the feature has no direction (strand==0), the returned patch will
        simply be a rectangle.

        The x-coordinates of the patch are determined by the feature's
        `start` and `end` while the y-coordinates are determined by the `level`.
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
        if is_undirected or head_is_cut:
            head_length = 0.001
        else:
            width_pixel = self._get_ax_width(ax, unit="pixel")
            head_length = 0.5 * width_pixel * feature.length / self.sequence_length
            head_length = min(head_length, 0.6 * feature.thickness)

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
        return "black" if (luminosity >= 0.5) else "white"

    def annotate_feature(
        self,
        ax,
        feature,
        level,
        inline=False,
        max_label_length=50,
        max_line_length=30,
        padding=0,
        indicate_strand_in_label=False,
    ):
        """Create a Matplotlib Text with the feature's label.

        The x-coordinates of the text are determined by the feature's
        `x_center` while the y-coordinates are determined by the `level`.

        The text is horizontally and vertically centered.

        The Arrow points in the direction of the feature's strand.
        If the feature has no direction (strand==0), the returned patch will
        simply be a rectangle.

        The x-coordinates of the patch are determined by the feature's
        `start` and `end` while the y-coordinates are determined by the `level`.
        """
        x, y = self.coordinates_in_plot(feature.x_center, level)
        label = feature.label
        if indicate_strand_in_label:
            if feature.strand == -1:
                label = "⇦" + label
            if feature.strand == 1:
                label = label + "⇨"

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
            bg_color = change_luminosity(feature.color, min_luminosity=0.95)
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

    def place_annotation(
        self,
        feature,
        ax,
        level,
        annotate_inline,
        max_line_length,
        max_label_length,
        indicate_strand_in_label=False,
    ):
        """"Place an annotation in the figure. Decide on inline vs. outline.

        Parameters
        ----------

        feature
          Graphic feature to place in the figure.

        ax
          Matplotlib ax in which to place the feature.

        level
          level at which the annotation should be placed.

        annotate_inline
          If true, the plotter will attempt to annotate inline, and fall back
          to outline annotation.

        max_line_length
          If an annotation label's length exceeds this number the label will
          wrap over several lines.

        max_label_length,
          If an annotation label's length exceeds this number the label will
          be cut with an ellipsis (...).

        indicate_strand_in_label
          If True, then the label will be represented as "<= label" or
          "label =>" with an arrow representing the strand.
        """
        padding = self.compute_padding(ax)
        if annotate_inline:
            # FIRST ATTEMPT TO ANNOTATE INSIDE THE FEATURE. CHECK FOR OVERFLOW
            text, overflowing, lines, (x1, x2), height = self.annotate_feature(
                ax=ax,
                feature=feature,
                level=level,
                inline=True,
                padding=padding,
                max_label_length=max_label_length,
                max_line_length=max_line_length,
                indicate_strand_in_label=indicate_strand_in_label,
            )

            # IF OVERFLOW, REMOVE THE TEXT AND PLACE IT AGAIN, OUTLINE.
            if overflowing:
                text.remove()
                text, _, lines, (x1, x2), height = self.annotate_feature(
                    ax=ax,
                    feature=feature,
                    level=level,
                    inline=False,
                    padding=padding,
                    max_label_length=max_label_length,
                    max_line_length=max_line_length,
                    indicate_strand_in_label=indicate_strand_in_label,
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
        strand_in_label_threshold="default",
        elevate_outline_annotations="default",
        x_lim=None,
        figure_height=None,
        sequence_params=None,
    ):
        """Plot all the features in the same Matplotlib ax.

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

        ruler_color
          Ruler color.

        plot_sequence
          If True and the graphic record has a "sequence" attribute set, the
          sequence will be displayed below the base line.

        annotate_inline
          If true, some feature labels will be displayed inside their
          corresponding feature if there is sufficient space.

        max_label_length
          If an annotation label's length exceeds this number the label will
          be cut with an ellipsis (...).

        max_line_length
          If an annotation label's length exceeds this number the label will
          wrap over several lines.

        level_offset
          All features and annotations will be pushed up by "level_offset". Can
          be useful when plotting several sets of features successively on a
          same ax.

        strand_in_label_pixel_threshold
          Number N such that, when provided, every feature with a graphical
          width in pixels below N will have its strand indicated in the label
          by an a left/right arrow.

        elevate_outline_annotations
          If true, every text annotation will be above every feature. If false,
          text annotations will be as close as possible to the features.

        x_lim
          Horizontal axis limits to be set at the end.

        figure_height
          Figure height.

        sequence_params
          parameters for plot_sequence.
        """

        if elevate_outline_annotations == "default":
            default = self.default_elevate_outline_annotations
            elevate_outline_annotations = default
        if strand_in_label_threshold == "default":
            default = self.default_strand_in_label_threshold
            strand_in_label_threshold = default

        features_levels = compute_features_levels(self.features)

        for f in features_levels:
            features_levels[f] += level_offset
        max_level = (
            1 if (features_levels == {}) else max(1, max(features_levels.values()))
        )
        auto_figure_height = (ax is None) and (figure_height is None)
        if ax is None:
            height = figure_height or max_level
            fig, ax = plt.subplots(1, figsize=(figure_width, height))

        def strand_in_label(f):
            """Anything under 0.1 inches in the figure."""
            if strand_in_label_threshold is None:
                return False
            width_pixel = self._get_ax_width(ax, unit="pixel")
            f_pixels = 1.0 * width_pixel * f.length / self.sequence_length
            return f_pixels < strand_in_label_threshold

        self.initialize_ax(ax, draw_line=draw_line, with_ruler=with_ruler)
        if x_lim is not None:
            ax.set_xlim(*x_lim)
        overflowing_annotations = []
        renderer = ax.figure.canvas.get_renderer()
        bbox = ax.get_window_extent(renderer)
        ax_height = bbox.height
        ideal_yspan = 0

        # sorting features from larger to smaller to make smaller features
        # appear "on top" of smaller ones, in case it happens. May be useless
        # now.
        sorted_features_levels = sorted(
            features_levels.items(), key=lambda o: -o[0].length
        )
        for feature, level in sorted_features_levels:
            self.plot_feature(ax=ax, feature=feature, level=level)
            if feature.label is None:
                continue
            (text, overflowing, nlines, (x1, x2,), height,) = self.place_annotation(
                feature=feature,
                ax=ax,
                level=level,
                annotate_inline=annotate_inline,
                max_line_length=max_line_length,
                max_label_length=max_label_length,
                indicate_strand_in_label=strand_in_label(feature),
            )
            line_height = height / nlines
            n_text_lines_in_axis = ax_height / line_height
            min_y_height = self.min_y_height_of_text_line
            feature_ideal_span = min_y_height * n_text_lines_in_axis
            ideal_yspan = max(ideal_yspan, feature_ideal_span)
            if overflowing or not annotate_inline:
                # trick here: we are representing text annotations as
                # GraphicFeatures so we can place them using
                # compute_features_levels().
                # We are also storing all the info necessary for label plotting
                # in these pseudo-graphic-features.
                overflowing_annotations.append(
                    GraphicFeature(
                        start=x1,
                        end=x2,
                        feature=feature,
                        text=text,
                        feature_level=level,
                        nlines=nlines,
                        color=feature.color,
                        label_link_color=feature.label_link_color,
                    )
                )

        # There are two ways to plot annotations: evelated, all above all the
        # graphic feature. Or at the same levels as the graphic features (
        # every annotation above its respective feature, but some annotations
        # can be below some features).
        if elevate_outline_annotations:

            base_feature = GraphicFeature(
                start=-self.sequence_length,
                end=self.sequence_length,
                fixed_level=0,
                nlines=1,
                is_base=True,
            )
            overflowing_annotations.append(base_feature)
            annotations_levels = compute_features_levels(overflowing_annotations)
        else:
            for f in self.features:
                f.data.update(dict(nlines=1, fixed_level=features_levels[f]))
            annotations_levels = compute_features_levels(
                overflowing_annotations + self.features
            )
            annotations_levels = {
                f: annotations_levels[f] for f in overflowing_annotations
            }

        max_annotations_level = max([0] + list(annotations_levels.values()))
        annotation_height = self.determine_annotation_height(max_level)
        annotation_height = max(self.min_y_height_of_text_line, annotation_height)
        labels_data = {}
        for feature, level in annotations_levels.items():
            if "is_base" in feature.data:
                continue
            text = feature.data["text"]
            x, y = text.get_position()
            if elevate_outline_annotations:
                new_y = (max_level) * self.feature_level_height + (
                    level
                ) * annotation_height
            else:
                new_y = annotation_height * level
            text.set_position((x, new_y))
            fx, fy = self.coordinates_in_plot(
                feature.data["feature"].x_center, feature.data["feature_level"]
            )

            # PLOT THE LABEL-TO-FEATURE LINK
            link_color = feature.label_link_color
            if link_color == "auto":
                link_color = change_luminosity(feature.color, luminosity=0.2)
            ax.plot([x, fx], [new_y, fy], c=link_color, lw=0.5, zorder=-10)
            labels_data[feature.data["feature"]] = dict(
                feature_y=fy, annotation_y=new_y
            )

        if plot_sequence:
            self.plot_sequence(ax, **(sequence_params or {}))

        self.finalize_ax(
            ax=ax,
            features_levels=max([1] + list(features_levels.values())),
            annotations_max_level=max_annotations_level,
            auto_figure_height=auto_figure_height,
            ideal_yspan=ideal_yspan,
            annotations_are_elevated=elevate_outline_annotations,
        )
        return ax, (features_levels, labels_data)

    def plot_legend(
        self, ax, allow_ambiguity=False, include_edge=True, **legend_kwargs
    ):
        handles = []
        features_parameters = {}
        for feature in self.features:
            text = feature.legend_text
            if text is None:
                continue
            parameters = dict(label=text, facecolor=feature.color, edgecolor="black",)
            if include_edge:
                parameters.update(
                    dict(linewidth=feature.linewidth, edgecolor=feature.linecolor,)
                )
            if text in features_parameters:
                previous_parameters = features_parameters[text]
                if (not allow_ambiguity) and any(
                    [parameters[k] != previous_parameters[k] for k in parameters]
                ):
                    raise ValueError("Cannot generate an unambiguous legend as two")
                continue
            features_parameters[text] = parameters
            handles.append(Patch(**parameters))
        ax.legend(handles=handles, **legend_kwargs)


def change_luminosity(color, luminosity=None, min_luminosity=None, factor=None):
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
      Only used if `luminosity` is not set. Positive factors increase
      luminosity and negative factors decrease it. More precisely, the
      luminosity of the new color is L^(-factor), where L is the current
      luminosity, between 0 and 1.
    """
    r, g, b = colorConverter.to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    new_l = l
    if luminosity is not None:
        new_l = luminosity
    if factor is not None:
        new_l = l ** (-factor)
    if min_luminosity is not None:
        new_l = max(new_l, min_luminosity)

    return colorsys.hls_to_rgb(h, new_l, s)


def get_text_box(text, margin=0):
    """Return the coordinates of a Matplotlib Text.

    `text` is a Matplotlib text obtained with ax.text().
    This returns `(x1,y1, x2, y2)` where (x1,y1) is the lower left corner
    and (x2, y2) is the upper right corner of the text, in data coordinates.
    If a margin m is supplied, the returned result is (x1-m, y1-m, x2+m, y2+m).
    """
    renderer = text.axes.figure.canvas.get_renderer()
    bbox = text.get_window_extent(renderer)  # bounding box
    __x1, y1, __x2, y2 = bbox.get_points().flatten()
    bbox = bbox.transformed(text.axes.transData.inverted())
    x1, __y1, x2, __y2 = bbox.get_points().flatten()
    return [x1, y1, x2, y2]
