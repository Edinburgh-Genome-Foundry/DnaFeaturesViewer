"""Simple example with plotly output. Requires the plotly library installed.
"""

from dna_features_viewer import BiopythonTranslator

record = BiopythonTranslator().translate_record(record="example_sequence.gb")
plot = record.plot_with_plotly(figure_width=8)

plot.write_html("plot_with_plotly.html", include_plotlyjs="cdn")
