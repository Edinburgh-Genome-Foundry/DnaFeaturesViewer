"""In this example we plot a record fragment with sequence over multiple lines.
"""
from dna_features_viewer import BiopythonTranslator


class CustomTranslator(BiopythonTranslator):
    def compute_feature_color(self, feature):
        return {
            "restriction_site": "yellow",
            "CDS": "orange",
            "promoter": "darkblue",
            "terminator": "lightblue",
        }[feature.type]


translator = CustomTranslator()
graphic_record = translator.translate_record("example_sequence.gb")
subrecord = graphic_record.crop((1800, 2750))
subrecord.plot_on_multiple_pages(
    "multipage_plot.pdf",
    nucl_per_line=70,
    lines_per_page=7,
    plot_sequence=True,
)
