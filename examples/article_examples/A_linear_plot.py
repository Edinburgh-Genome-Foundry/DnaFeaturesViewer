from dna_features_viewer import BiopythonTranslator


class CustomTranslator(BiopythonTranslator):

    # Label fields indicates the order in which annotations fields are
    # considered to determine the feature's label
    label_fields = ["label", "note", "name", "gene"]

    def compute_feature_legend_text(self, feature):
        return feature.type
    
    def compute_feature_color(self, feature):
        return {
            "rep_origin": "yellow",
            "CDS": "#ffd383",  # light orange
            "regulatory": "red",
            "misc_recomb": "#fbf3f6",  # pink
            "misc_feature": "#d1e9f1",  # light blue
            "backbone": "darkblue",
        }[feature.type]
    
    def compute_feature_box_color(self, feature):
        return "white"
    
    def compute_feature_box_linewidth(self, feature):
        return 0


translator = CustomTranslator()
graphic_record = translator.translate_record("plasmid.gb")
ax, _ = graphic_record.plot(figure_width=13, strand_in_label_threshold=7)
graphic_record.plot_legend(ax=ax, loc=1, ncol=3, frameon=False)
ax.figure.savefig("A_linear_plot.svg", bbox_inches="tight")
