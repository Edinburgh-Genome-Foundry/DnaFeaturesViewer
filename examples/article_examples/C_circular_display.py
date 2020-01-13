from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord

class ExpressionUnitTranslator(BiopythonTranslator):

    def compute_feature_color(self, feature):
        color_map = {
            "rep_origin": "yellow",
            "CDS": "#ffd383",  # light orange
            "regulatory": "red",
            "misc_recomb": "darkblue",
            "misc_feature": "#d1e9f1",  # light blue
            "backbone": "darkblue",
        }
        return color_map[feature.type]

    def compute_feature_label(self, feature):
        if feature.type not in ["CDS", "regulatory"]:
            return None
        else:
            return BiopythonTranslator.compute_feature_label(self, feature)


translator = ExpressionUnitTranslator()
graphic_record = translator.translate_record(
    "plasmid.gb", record_class=CircularGraphicRecord
)
graphic_record.top_position = 4800  # sequence index appearing at the top
ax, _ = graphic_record.plot(figure_width=4)
ax.figure.savefig("C_circular_display.svg", bbox_inches="tight")