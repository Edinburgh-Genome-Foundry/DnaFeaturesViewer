try:
    import plotly.graph_objects as go

    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

try:
    import pandas as pd

    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

import matplotlib.pyplot as plt


class PlotlyPlottableMixin:
    def plotly_feature_patch(
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
        """Return a dict with points coordinates of a plotly shape. Same as bokeh feature arrow

        Parameters
        ----------

        start, end, strand

        """
        hw = width / 2.0
        x1, x2 = (start, end) if (strand >= 0) else (end, start)
        bp_per_width = figure_width / self.sequence_length
        delta = arrow_width_inches / bp_per_width
        if strand >= 0:
            head_base = max(x1, x2 - delta)
        else:
            head_base = min(x1, x2 + delta)
        result = dict(
            xs=[x1, x1, head_base, x2, head_base, x1],
            ys=[e + level for e in [-hw, hw, hw, 0, -hw, -hw]],
        )
        result.update(kwargs)
        return result

    def plot_with_plotly(self, figure_width=5, figure_height="auto", tools="auto"):
        """Plot the graphic record using Plotly.
        The returned fig object can be used in a Dash dashboard.

        Examples
        --------

        >>>


        """
        if not PLOTLY_AVAILABLE:
            raise ImportError("``plot_with_plotly`` requires plotly installed.")
        if not PANDAS_AVAILABLE:
            raise ImportError("``plot_with_plotly`` requires plotly installed.")

        # Set up default tools
        # if tools == "auto":
        #     tools = [HoverTool(tooltips="@hover_html"), "xpan,xwheel_zoom,reset,tap"]

        # FIRST PLOT WITH MATPLOTLIB AND GATHER INFOS ON THE PLOT
        ax, (features_levels, plot_data) = self.plot(figure_width=figure_width)
        width, height = [int(100 * e) for e in ax.figure.get_size_inches()]
        plt.close(ax.figure)
        if figure_height == "auto":
            height = int(0.5 * height)
        else:
            height = 100 * figure_height
        height = max(height, 185) # Minimal height to see all icons

        max_y = max(
            [data["annotation_y"] for f, data in plot_data.items()]
            + list(features_levels.values())
        )

        fig = go.Figure()

        # Update plot width and height
        fig.update_layout(
            autosize=False,
            width=width,
            height=height,
            margin=dict(l=0, r=20, t=0, b=20)
        )

        # Update axes properties
        fig.update_xaxes(
            range=[0, self.sequence_length],
            zeroline=False,
        )

        fig.update_yaxes(
            range=[-1, max_y + 1],
            zeroline=False,
            showline=False,
            showgrid=False,
            visible=False,
        )

        # Add patches
        for feature, level in features_levels.items():
            patch = self.plotly_feature_patch(
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
            fig.add_trace(
                go.Scatter(
                    x=patch["xs"],
                    y=patch["ys"],
                    fill="toself",
                    mode="lines",
                    name="",
                    text=patch["label"],
                    line_color="#000000",
                    fillcolor=patch["color"],
                )
            )

        if plot_data != {}:

            text_df = pd.DataFrame.from_records(
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

            # Scatter trace of text labels
            fig.add_trace(go.Scatter(
                x=text_df["x"],
                y=text_df["y"],
                text=text_df["text"],
                mode="text",
                hoverinfo="skip",
                textfont=dict(size=12, family="arial"),
                textposition="middle center",
            ))

            # Add segments
            for feature, pdata in plot_data.items():
                fig.add_shape(
                    type="line",
                    x0=feature.x_center,
                    y0=pdata["annotation_y"],
                    x1=feature.x_center,
                    y1=pdata["feature_y"],
                    line=dict(color="#000000", width=0.5)
                )

        fig.update_layout(
            template="simple_white",
            showlegend=False,
        )

        return fig
