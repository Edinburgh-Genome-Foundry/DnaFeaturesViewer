import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

from .tools import (change_luminosity, get_text_box, compute_features_levels,
                    bokeh_feature_patch)
from .biotools import extract_translation
from .GraphicFeature import GraphicFeature

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Alphabet import DNAAlphabet


try:
    from bokeh.plotting import figure, ColumnDataSource
    from bokeh.models import Range1d, HoverTool
    BOKEH_AVAILABLE = True
except:
    BOKEH_AVAILABLE = False

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except:
    PANDAS_AVAILABLE = False


class GraphicRecord:
    """Set of Genetic Features of a same DNA sequence, to be plotted together.

    Parameters
    ----------

    sequence_length
      Length of the DNA sequence, in number of nucleotides

    features
      list of GraphicalFeature objects.

    feature_level_height
      Width in inches of one "level" for feature arrows.

    annotation_height
      Width in inches of one "level" for feature annotations.

    first_index
      Indicates the first index to plot in case the sequence is actually a
      subsequence of a larger one. For instance, if the Graphic record
      represents the segment (400, 420) of a sequence, we will have
      ``first_index=400`` and ``sequence_length=20``.

    plots_indexing
      Indicates which standard to use to show nucleotide indices in the plots.
      If 'biopython', the standard python indexing is used (starting at 0).
      If 'genbank', the indexing follows the Genbank standard (starting at 1).
    """

    def __init__(self, sequence_length=None, sequence=None, features=(),
                 feature_level_height=1,
                 first_index=0, plots_indexing='biopython',
                 labels_spacing=10):
        if sequence_length is None:
            sequence_length = len(sequence)
        self.features = features
        self.sequence_length = sequence_length
        self.feature_level_height = feature_level_height
        self.sequence = sequence
        self.first_index = first_index
        self.plots_indexing = plots_indexing
        self.labels_spacing = labels_spacing

    @property
    def span(self):
        """Return the display span (start, end) accounting for first_index."""
        return self.first_index, self.first_index + self.sequence_length

    def initialize_ax(self, ax, draw_line, with_ruler):
        """Initialize the ax: remove axis, draw a horizontal line, etc.

        Parameters
        ----------
        
        draw_line
          True/False to draw the horizontal line or not.

        with_ruler
          True/False to draw the indices indicators along the line. 

        """
        start, end = self.span
        plot_start, plot_end = start - 0.8, end - 0.2
        if draw_line:
            ax.plot([plot_start, plot_end], [0, 0], zorder=-1000, c="k")

        if with_ruler:  # only display the xaxis ticks
            ax.set_frame_on(False)
            ax.yaxis.set_visible(False)
            ax.xaxis.tick_bottom()
        else:  # don't display anything
            ax.axis("off")

        ax.set_xlim(plot_start, plot_end)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

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
        head_is_cut = ((feature.strand == 1 and feature.open_right) or
                       (feature.strand == -1 and feature.open_left))
        head_length = (0.001 if (is_undirected or head_is_cut) else
                       max(0.6 * feature.thickness, 5))

        arrowstyle = mpatches.ArrowStyle.Simple(head_width=feature.thickness,
                                                tail_width=feature.thickness,
                                                head_length=head_length)
        y = self.feature_level_height * level
        patch = mpatches.FancyArrowPatch([x1, y], [x2, y],
                                         shrinkA=0.0, shrinkB=0.0,
                                         arrowstyle=arrowstyle,
                                         facecolor=feature.color, zorder=0,
                                         edgecolor=feature.linecolor,
                                         linewidth=feature.linewidth)
        ax.add_patch(patch)
        return patch

    def coordinates_in_plot(self, x, level):
        """Convert a sequence position and height level into a (x, y) position.
        """
        return (x, level * self.feature_level_height)

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
                f1, f2 = f.split_in_two(self.sequence_length - 1)
                f2.start, f2.end = (f2.start - self.sequence_length,
                                    f2.end - self.sequence_length)
                new_features += [f1, f2]
            else:
                new_features.append(f)
        self.features = new_features

    def annotate_feature(self, ax, feature, level):
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
        box_color = feature.box_color
        text = ax.text(
            x, y, feature.label,
            horizontalalignment="center",
            verticalalignment="center",
            bbox=None if (box_color is None) else dict(boxstyle="round",
                fc=bg_color if box_color is "auto" else box_color,
                ec="0.5", lw=feature.box_linewidth),
            fontdict=feature.fontdict,
            zorder=2
        )
        figure_width = ax.figure.get_size_inches()[0]
        x1, y1, x2, y2 = get_text_box(text)
        x1 -= self.labels_spacing
        x2 += self.labels_spacing
        overflowing = (x1 < feature.start) or (x2 > feature.end)
        return text, overflowing, (x1, x2), (y2 - y1)
    
    def determine_annotation_height(self, levels):
        """By default the ideal annotation level height is the same as the
        feature_level_height."""
        # Improve me: ideally, annotation width would be linked to the height
        # of one line of text, so dependent on font size and ax height/span.  
        return self.feature_level_height

    def plot(self, ax=None, figure_width=8, draw_line=True, with_ruler=True,
             plot_sequence=False, annotate_inline=False, level_offset=0,
             x_lim=None):
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
        max_level = (1 if (features_levels == {}) else
                     max(1, max(features_levels.values())))
        auto_figure_height = ax is None
        if ax is None:
            fig, ax = plt.subplots(1, figsize=(figure_width, max_level))
        self.initialize_ax(ax, draw_line=draw_line, with_ruler=with_ruler)
        if x_lim is not None:
            ax.set_xlim(*x_lim)
        bbox = ax.get_window_extent()
        ax_width, ax_height = bbox.width, bbox.height

        overflowing_annotations = []
        ideal_yspan = 0
        for feature, level in features_levels.items():
            self.plot_feature(ax=ax, feature=feature, level=level)
            if feature.label is not None:
                text, overflowing, (x1, x2), height = self.annotate_feature(
                    ax=ax, feature=feature, level=level
                )
                nlines = len(feature.label.split("\n"))
                line_height = height / nlines
                # Hey look, a magic number!
                feature_ideal_span = 0.4 * (ax_height / line_height)
                ideal_yspan = max(ideal_yspan,  feature_ideal_span)
                if overflowing or not annotate_inline:
                    overflowing_annotations.append(GraphicFeature(
                        start=x1, end=x2, feature=feature,
                        text=text, feature_level=level,
                        nlines=nlines
                    ))
        annotations_levels = compute_features_levels(overflowing_annotations)
        max_annotations_level = max([0] + list(annotations_levels.values()))
        annotation_height = self.determine_annotation_height(max_level)
        labels_data = {}
        for feature, level in annotations_levels.items():
            text = feature.data["text"]
            x, y = text.get_position()
            new_y = ((max_level + 1) * self.feature_level_height +
                     (level) * annotation_height)
            text.set_position((x, new_y))
            fx, fy = self.coordinates_in_plot(feature.data["feature"].x_center,
                                              feature.data["feature_level"])
            ax.plot([x, fx], [new_y, fy], c="k", lw=0.5, zorder=1)
            labels_data[feature.data["feature"]] = dict(
                feature_y=fy,
                annotation_y=new_y
            )

        if plot_sequence:
            self.plot_sequence(ax)
        
        self.finalize_ax(
            ax=ax,
            features_levels=max([1] + list(features_levels.values())),
            annotations_max_level=max_annotations_level,
            auto_figure_height=auto_figure_height,
                         ideal_yspan=ideal_yspan)
        return ax, (features_levels, labels_data)

    def plot_sequence(self, ax, location=None, y_offset=1, fontdict=None,
                      background=("#f7fbff", "#fffcf0")):
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
        lstart, lend = location
        fontdict = dict(size=11)
        fontdict.update(fontdict or {})
        for i, n in enumerate(self.sequence):
            l = i + lstart
            if (lstart <= l <= lend):
                ax.text(l, - 0.7 * self.feature_level_height * y_offset,
                        n, ha='center', va='center', fontdict=fontdict)
        if background is not None:
            for i in range(lstart, lend):
                ax.fill_between([i - 0.5, i + 0.5], y1=1000, y2=-1000,
                                zorder=-2000, facecolor=background[i % 2])
        ymin = ax.get_ylim()[0]
        ax.set_ylim(bottom=min(ymin, -y_offset * self.feature_level_height))

    def plot_translation(self, ax, location=None, y_offset=2, fontdict=None,
                         background=("#f5fff0", "#fff7fd"), translation=None,
                         long_form_translation=True):
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
            translation = extract_translation(self.sequence, location=new_loc,
                                              long_form=long_form_translation)
        texts = [
            ((start + 3 * i, start + 3 * (i + 1)), aa)
            for i, aa in enumerate(translation)
        ]

        y = - 0.7 * y_offset * self.feature_level_height
        ymin = ax.get_ylim()[0]
        ax.set_ylim(bottom=min(ymin, -y_offset * self.feature_level_height))
        fontdict = fontdict or {}
        for i, ((start, end), text) in enumerate(texts):
            ax.text(0.5 * (start + end - 1), y, text,
                    ha='center', va='center', fontdict=fontdict)
            if background:
                ax.fill_between([start - 0.5, end - 0.5], y1=1000, y2=-1000,
                                zorder=-1000, facecolor=background[i % 2])

    def finalize_ax(self, ax, features_levels, annotations_max_level,
                    auto_figure_height=False,
                    ideal_yspan=None):
        annotation_height = self.determine_annotation_height(None)
        ymax = (self.feature_level_height * (features_levels + 1) +
                annotation_height * (annotations_max_level + 1))
        ymin = -1
        if ideal_yspan is not None:
            ymax = max(ideal_yspan + ymin, ymax)
        ax.set_ylim(ymin, ymax)
        if auto_figure_height:
            figure_width = ax.figure.get_size_inches()[0]
            ax.figure.set_size_inches(figure_width, 1 + 0.4 * ymax)
        if self.plots_indexing == 'genbank':
            ax.set_xticklabels([int(i + 1) for i in ax.get_xticks()])
        return ideal_yspan / (ymax - ymin)

    def to_biopython_record(self, sequence):
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
        features = [
            SeqFeature(FeatureLocation(f.start, f.end, f.strand),
                       type=f.feature_type, qualifiers={"label": f.label})
            for f in self.features
        ]
        if not isinstance(sequence, Seq):
            sequence = Seq(sequence, alphabet=DNAAlphabet())
        return SeqRecord(seq=sequence, features=features)

    def crop(self, window):
        s, e = window
        if (s < 0) or (e >= self.sequence_length):
            raise ValueError("out-of-bound cropping")
        new_features = []
        for f in self.features:
            cropped_feature = f.crop(window)
            if cropped_feature is not None:  # = has ovelap with the window
                new_features.append(cropped_feature)

        return GraphicRecord(
            sequence=self.sequence[s:e] if self.sequence is not None else None,
            sequence_length=e - s,
            features=new_features,
            feature_level_height=self.feature_level_height,
            first_index=self.first_index + s
        )

    def plot_with_bokeh(self, figure_width=5):
        """Plot the graphic record using Bokeh.

        Examples
        --------

        >>>


        """
        if not BOKEH_AVAILABLE:
            raise ImportError("``plot_with_bokeh`` requires Bokeh installed.")
        if not PANDAS_AVAILABLE:
            raise ImportError("``plot_with_bokeh`` requires Pandas installed.")
        ax, (features_levels, plot_data) = self.plot(figure_width=figure_width)
        width, height = [int(100 * e) for e in ax.figure.get_size_inches()]
        plt.close(ax.figure)
        max_y = max([data["annotation_y"] for f, data in plot_data.items()] +
                    list(features_levels.values()))
        hover = HoverTool(tooltips="@hover_html")
        p = figure(plot_width=width, plot_height=height,
                   tools=[hover, "xpan,xwheel_zoom,reset,tap"],
                   x_range=Range1d(0, self.sequence_length),
                   y_range=Range1d(-1, max_y + 1))
        p.patches(
            xs='xs', ys='ys', color='color', line_color="#000000",
            source=ColumnDataSource(pd.DataFrame.from_records([
                bokeh_feature_patch(
                    self, feature.start, feature.end, feature.strand,
                    level=level, color=feature.color,
                    label=feature.label,
                    hover_html=(feature.html if feature.html is not None else
                                feature.label)
                )
                for feature, level in features_levels.items()
            ]))
        )

        if plot_data != {}:
            p.text(
                x='x', y='y', text='text', text_align="center",
                text_font_size='12px',  text_font="arial",
                text_font_style="normal",
                source=ColumnDataSource(pd.DataFrame.from_records([
                    dict(x=feature.x_center, y=pdata["annotation_y"],
                         text=feature.label, color=feature.color)
                    for feature, pdata in plot_data.items()
                ]))
            )
            p.segment(
                x0='x0', x1='x1', y0='y0', y1='y1', line_width=0.5,
                color="#000000",
                source=ColumnDataSource(pd.DataFrame.from_records([
                    dict(x0=feature.x_center, x1=feature.x_center,
                         y0=pdata["annotation_y"], y1=pdata["feature_y"])
                    for feature, pdata in plot_data.items()
                ]))
            )

        p.yaxis.visible = False
        p.outline_line_color = None
        p.grid.grid_line_color = None
        p.toolbar.logo = None

        return p
