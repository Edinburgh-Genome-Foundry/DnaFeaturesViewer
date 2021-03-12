try:
    import bokeh
    from bokeh.plotting import figure, ColumnDataSource
    from bokeh.models import Range1d, HoverTool
    from bokeh.core.properties import value

    BOKEH_AVAILABLE = True
except ImportError:
    BOKEH_AVAILABLE = False

try:
    import pandas as pd

    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

from packaging import version
import matplotlib.pyplot as plt


class BokehPlottableMixin:
    def bokeh_feature_patch(
        self,
        start,
        end,
        strand,
        figure_width=5,
        width=0.4,
        level=0,
        arrow_width_inches=0.05,
        **kwargs
    ):
        """Return a dict with points coordinates of a Bokeh Feature arrow.

        Parameters
        ----------

        start, end, strand

        """
        hw = width / 2.0
        x1, x2 = (start, end) if (strand >= 0) else (end, start)
        bp_per_width = figure_width / self.sequence_length
        delta = arrow_width_inches / bp_per_width
        if strand > 0:
            head_base = max(x1, x2 - delta)
        elif strand < 0:
            head_base = min(x1, x2 + delta)
        else:
            head_base = x2
        result = dict(
            xs=[x1, x1, head_base, x2, head_base, x1],
            ys=[e + level for e in [-hw, hw, hw, 0, -hw, -hw]],
        )
        result.update(kwargs)
        return result

    def plot_with_bokeh(self, figure_width=5, figure_height="auto", tools="auto"):
        """Plot the graphic record using Bokeh.

        Examples
        --------

        >>>


        """
        if not BOKEH_AVAILABLE:
            raise ImportError("``plot_with_bokeh`` requires Bokeh installed.")
        if not PANDAS_AVAILABLE:
            raise ImportError("``plot_with_bokeh`` requires Pandas installed.")

        # Set up default tools
        if tools == "auto":
            tools = [HoverTool(tooltips="@hover_html"), "xpan,xwheel_zoom,reset,tap"]

        # FIRST PLOT WITH MATPLOTLIB AND GATHER INFOS ON THE PLOT
        ax, (features_levels, plot_data) = self.plot(figure_width=figure_width)
        width, height = [int(100 * e) for e in ax.figure.get_size_inches()]
        plt.close(ax.figure)
        if figure_height == "auto":
            height = int(0.5 * height)
        else:
            height = 100 * figure_height
        height = max(height, 185)  # Minimal height to see all icons

        max_y = max(
            [data["annotation_y"] for f, data in plot_data.items()]
            + list(features_levels.values())
        )

        # BUILD THE PLOT ()
        plot = figure(
            plot_width=width,
            plot_height=height,
            tools=tools,
            x_range=Range1d(0, self.sequence_length),
            y_range=Range1d(-1, max_y + 1),
        )
        plot.patches(
            xs="xs",
            ys="ys",
            color="color",
            line_color="#000000",
            source=ColumnDataSource(
                pd.DataFrame.from_records(
                    [
                        self.bokeh_feature_patch(
                            feature.start,
                            feature.end,
                            feature.strand,
                            figure_width=figure_width,
                            level=level,
                            color=feature.color,
                            label=feature.label,
                            hover_html=(
                                feature.html
                                if feature.html is not None
                                else feature.label
                            ),
                        )
                        for feature, level in features_levels.items()
                    ]
                )
            ),
        )

        if plot_data != {}:
            if version.parse(bokeh.__version__) < version.parse("2.3"):
                value_arial = "arial"
            else:  # >= 2.3
                value_arial = value("arial")
            plot.text(
                x="x",
                y="y",
                text="text",
                text_align="center",
                text_font_size="12px",
                text_font=value_arial,
                text_font_style="normal",
                source=ColumnDataSource(
                    pd.DataFrame.from_records(
                        [
                            dict(
                                x=feature.x_center,
                                y=pdata["annotation_y"],
                                text=feature.label,
                                color=feature.color,
                            )
                            for feature, pdata in plot_data.items()
                        ]
                    )
                ),
            )
            plot.segment(
                x0="x0",
                x1="x1",
                y0="y0",
                y1="y1",
                line_width=0.5,
                color="#000000",
                source=ColumnDataSource(
                    pd.DataFrame.from_records(
                        [
                            dict(
                                x0=feature.x_center,
                                x1=feature.x_center,
                                y0=pdata["annotation_y"],
                                y1=pdata["feature_y"],
                            )
                            for feature, pdata in plot_data.items()
                        ]
                    )
                ),
            )

        plot.yaxis.visible = False
        plot.outline_line_color = None
        plot.grid.grid_line_color = None
        plot.toolbar.logo = None

        return plot
