"""Implements the missing Matplotlib ArrowWedge patch class.

This is a plain arrow curved alongside a protion of circle, like you would
expect a circular genetic feature to look.
"""
import numpy as np
import matplotlib.patches as mpatches


class ArrowWedge(mpatches.Wedge):
    """Matplotlib patch shaped as a tick fraction of circle with a pointy end.

    This is the patch used by CircularGraphicRecord to draw features.

    Parameters
    ----------

    center
      Center of the circle around which the arrow-wedge is drawn.

    radius
      Radius of the circle around which the arrow-wedge is drawn.

    theta1
      Start angle of the wedge.

    theta2
      End angle of the wedge.

    width
      Width or thickness of the arrow-wedge.

    direction
      Determines whether the pointy end points in direct sense (+1) or
      indirect sense (-1) or no sense at all (0).
    """

    def __init__(self, center, radius, theta1, theta2, width, direction=+1, **kwargs):
        self.direction = direction
        self.radius = radius
        mpatches.Wedge.__init__(
            self, center, radius, theta1, theta2, width=width, **kwargs
        )
        self._recompute_path()

    def _recompute_path(self):
        """Recompute the full path forming the "tick" arrowed wedge.

        This method overwrites "mpatches.Wedge._recompute_path" in the
        super-class.
        """

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
                np.array([np.cos(np.deg2rad(theta1)), np.sin(np.deg2rad(theta1))]),
                inner_arc[0],
            ]
        else:
            angle_start_arrow = theta2 - arrow_angle
            arc = mpatches.Path.arc(theta1, angle_start_arrow)
            outer_arc = arc.vertices * (self.radius + self.width / 2.0) / self.radius
            inner_arc = (
                arc.vertices[::-1] * (self.radius - self.width / 2.0) / self.radius
            )
            arrow_vertices = [
                outer_arc[-1],
                np.array([np.cos(np.deg2rad(theta2)), np.sin(np.deg2rad(theta2))]),
                inner_arc[0],
            ]
        p = np.vstack([outer_arc, arrow_vertices, inner_arc])

        path_vertices = np.vstack([p, inner_arc[-1, :], (0, 0)])

        path_codes = np.hstack(
            [
                arc.codes,
                4 * [mpatches.Path.LINETO],
                arc.codes[1:],
                mpatches.Path.LINETO,
                mpatches.Path.CLOSEPOLY,
            ]
        )
        path_codes[len(arc.codes)] = mpatches.Path.LINETO

        # Shift and scale the wedge to the final location.
        path_vertices *= self.r
        path_vertices += np.asarray(self.center)
        self._path = mpatches.Path(path_vertices, path_codes)
