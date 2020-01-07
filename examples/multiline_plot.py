"""In this example we plot a record fragment with sequence over multiple lines.
"""
from dna_features_viewer import BiopythonTranslator

translator = BiopythonTranslator()
graphic_record = translator.translate_record("example_sequence.gb")
subrecord = graphic_record.crop((1700, 2000))
fig, axes = subrecord.plot_on_multiple_lines(
    nucl_per_line=70, plot_sequence=True
)
fig.savefig("multiline_plot.png")
