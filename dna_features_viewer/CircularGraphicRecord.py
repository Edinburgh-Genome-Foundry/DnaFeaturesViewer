import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from .GraphicRecord import GraphicRecord
import numpy as np

class ArrowWedge(mpatches.Wedge):

    def __init__(self, center, radius, theta1, theta2, width, direction=+1,
                 **kwargs):

        self.direction = direction
        self.radius = radius
        mpatches.Wedge.__init__(self, center, radius,
                                theta1, theta2, width,
                                **kwargs)
        self._recompute_path()

    def _recompute_path(self):

        if self.direction not in [-1, +1]:
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
    """Set of Genetic Features of a same DNA sequence, to be plotted together.

    Parameters
    ----------

    sequence_length
      Length of the DNA sequence, in number of nucleotides

    features
      list of GraphicalFeature objects.

    top_position
      The index in the sequence that will end up at the top of the circle

    feature_level_width
      Width in inches of one "level" for feature arrows.

    annotation_level_width
      Width in inches of one "level" for feature annotations.
    """

    def __init__(self, sequence_length, features, top_position=0,
                 feature_level_width=0.2, annotation_level_width=0.25, **kw):

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
                (0, -self.radius), self.radius, facecolor='none',
                edgecolor='k')
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
                           edgecolor='k',
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
