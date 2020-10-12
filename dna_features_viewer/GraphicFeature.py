from copy import deepcopy


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
      Short descriptive text associated and plotted with the feature.

    locus_tag
      Locus tag of the feature

    color
      Color of the feature, any Matplotlib-compatible format is accepted,
      such as "white", "w", "#ffffff", (1,1,1), etc.

    linecolor
      Color of the feature's border, any Matplotlib-compatible format is
      accepted, such as "white", "w", "#ffffff", (1,1,1), etc.

    linewidth
      Width of the line (=edge) surrounding the graphic feature, in pixels.

    thickness
      Vertical span of the feature.

    box_color
      Color of the label box. Set to None for no box around the label.
      Leave to "auto" for a box color that is a lightened version of the
      feature's color.

    data
      Any other keyword is kept into the feature.data[] dictionary.

    fontdict
      A Matplotlib fontdict for the font to be used in the label, e.g.
      ``size=11``, ``weight='bold'``, ``family='Helvetica'``, etc.

    open_left, open_right
      Set to True if this feature does not end on the right or left because it
      is a cropped version of a bigger feature.

    box_linewidth
      Width of the line delimiting the text box when the annotation is outside
      the graphic feature. Set to 0 for no box borders.

    box_color
      Background color of the annotation's text box. If left to "auto" the
      color will be a lighter version of the feature's color.

    label_link_color
      Color of the line linking the text annotation to its respective graphic
      feature. Set to auto for the line to automatically be a darker version
      of the feature's color.
    """

    feature_type = "feature"

    def __init__(
        self,
        start=None,
        end=None,
        strand=None,
        label=None,
        locus_tag=None,
        color="#000080",
        thickness=14,
        linewidth=1.0,
        linecolor="#000000",
        fontdict=None,
        html=None,
        open_left=False,
        open_right=False,
        box_linewidth=1,
        box_color="auto",
        legend_text=None,
        label_link_color="black",
        **data
    ):
        self.start = start
        self.end = end
        self.strand = strand
        self.label = label
        self.locus_tag = locus_tag
        self.color = color
        self.linecolor = linecolor
        self.data = data
        self.thickness = thickness
        self.linewidth = linewidth
        self.box_linewidth = box_linewidth
        self.box_color = box_color
        self.label_link_color = label_link_color
        self.fontdict = dict([("fontsize", 11)] + list((fontdict or {}).items()))
        self.html = html
        self.open_left = open_left
        self.open_right = open_right
        self.legend_text = legend_text

    def split_in_two(self, x_coord=0):
        """Return two features by cutting this feature at x_coord."""
        copy1 = deepcopy(self)
        copy2 = deepcopy(self)
        copy1.end = x_coord
        copy2.start = x_coord + 1
        return copy1, copy2

    def crop(self, window):
        """Return a the fragment of the feature that is in the window.

        If there is no overlap between the feature location and the window,
        None is returned.
        """
        s, e = window
        if (s > self.end) or (e < self.start):
            return None
        copy = deepcopy(self)
        if s > self.start:
            copy.start = s
            copy.open_left = True
        if e < self.end:
            copy.end = e
            copy.open_right = True
        return copy

    def overlaps_with(self, other):
        """Return whether the feature's location overlaps with feature `other`.
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
        return 0.5 * (self.start + self.end - 1)

    @staticmethod
    def from_biopython_feature(feature, **props):
        """Create a GraphicalFeature from a Biopython.Feature object."""
        return GraphicFeature(
            start=feature.location.start,
            end=feature.location.end,
            strand=feature.location.strand,
            **props
        )

    def __repr__(self):
        return ("GF(%(label)s, %(start)d-%(end)d " % self.__dict__) + (
            ")" if self.strand is None else "(%d))" % self.strand
        )
