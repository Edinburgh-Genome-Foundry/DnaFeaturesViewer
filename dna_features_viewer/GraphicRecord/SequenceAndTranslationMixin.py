from ..biotools import extract_graphical_translation


class SequenceAndTranslationMixin:
    def plot_sequence(
        self, ax, location=None, y_offset=1, fontdict=None, guides_intensity=0
    ):
        """Plot a sequence of nucleotides at the bottom of the plot.

        Parameters
        ----------

        ax
          Which axes the translation should be plotted to.

        location
          location of the segment to translate, either (start, end) or
          (start, end, strand).

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

        guides_intensity
          Intensity of the vertical guides marking the different nucleotides
          (0 = no guides).
        """
        if self.sequence is None:
            raise ValueError("No sequence in the graphic record")
        if location is None:
            location = self.span
        location_start, location_end = location
        fontdict = {"size": 11, **(fontdict or {})}
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
        if guides_intensity:
            color = (0, 0, 0, guides_intensity)
            for i in range(location_start, location_end + 1):
                ax.axvline(i - 0.5, linewidth=0.1, color=color, zorder=-10000)
        ymin = ax.get_ylim()[0]
        if ymin < -500:
            ymin = 0
        ax.set_ylim(bottom=min(ymin, -y_offset * self.feature_level_height))

    def plot_translation(
        self,
        ax,
        location=None,
        y_offset=2,
        fontdict=None,
        guides_intensity=0.5,
        translation=None,
        long_form_translation=True,
    ):
        """Plot a sequence of amino-acids at the bottom of the plot.

        Parameters
        ----------

        ax
          Which axes the translation should be plotted to.

        location
          location of the segment to translate (start, end).

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
                self.sequence, location=new_loc, long_form=long_form_translation,
            )
        texts = [
            ((start + 3 * i, start + 3 * (i + 1)), aa)
            for i, aa in enumerate(translation)
        ]

        y = -0.7 * y_offset * self.feature_level_height
        ymin = ax.get_ylim()[0]
        ax.set_ylim(bottom=min(ymin, -y_offset * self.feature_level_height))
        fontdict = {"size": 11, **(fontdict or {})}
        guides_color = (0, 0, 0, guides_intensity)
        for i, ((start, end), text) in enumerate(texts):
            ax.text(
                0.5 * (start + end - 1),
                y,
                text,
                ha="center",
                va="center",
                fontdict=fontdict,
            )
            if guides_intensity:
                ax.axvline(
                    start - 0.5, linewidth=0.1, color=guides_color, zorder=-10000,
                )
        if guides_intensity:
            ax.axvline(end - 0.5, linewidth=0.1, color=guides_color, zorder=-10000)
