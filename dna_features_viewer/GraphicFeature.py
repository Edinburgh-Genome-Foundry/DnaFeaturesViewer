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
                 html=None, open_left=False, open_right=False, **data):
        self.start = start
        self.end = end
        self.strand = strand
        self.label = label
        self.color = color
        self.data = data
        self.thickness = thickness
        self.linewidth = linewidth
        self.html = html
        self.open_left = open_left
        self.open_right = open_right

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
