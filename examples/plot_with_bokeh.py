"""Simple example with Bokeh output. Requires the Bokeh library installed.
"""

from dna_features_viewer import BiopythonTranslator
from bokeh.resources import CDN
from bokeh.embed import file_html

record = BiopythonTranslator().translate_record(record="example_sequence.gb")
plot = record.plot_with_bokeh(figure_width=8)

with open("plot_with_bokeh.html", "w+") as f:
    f.write(file_html(plot, CDN, "Example Sequence"))
