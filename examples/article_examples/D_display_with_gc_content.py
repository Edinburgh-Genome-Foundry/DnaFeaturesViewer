from dna_features_viewer import BlackBoxlessLabelTranslator, load_record
import matplotlib.pyplot as plt
import numpy as np


class GCIndicatingTranslator(BlackBoxlessLabelTranslator):
    def compute_feature_legend_text(self, feature):
        if feature.qualifiers["gc%"] < 30:
            return "GC < 30%"
        elif feature.qualifiers["gc%"] < 60:
            return "30-60% GC"
        else:
            return "GC > 60%"

    def compute_feature_color(self, feature):
        return {
            "GC < 30%": "peachpuff",
            "30-60% GC": "azure",
            "GC > 60%": "skyblue",
        }[self.compute_feature_legend_text(feature)]

    def compute_feature_fontdict(self, feature):
        return dict(size=10, weight="bold", color="#494949")

    def compute_feature_label(self, feature):
        if not (30 < feature.qualifiers["gc%"] < 60):
            normal_label = super().compute_feature_label(feature)
            return normal_label + "-%d%%" % feature.qualifiers["gc%"]


def gc_content(sequence):
    return 100.0 * len([c for c in sequence if c in "GC"]) / len(sequence)


# DISPLAY THE SEQUENCE MAP

fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(6, 3.5), sharex=True, gridspec_kw={"height_ratios": [3, 1]},
)
record = load_record("plasmid.gb")
for feature in record.features:
    feature.qualifiers["gc%"] = gc_content(feature.location.extract(record))
translator = GCIndicatingTranslator()
graphic_record = translator.translate_record(record)

graphic_record.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=7)
graphic_record.plot_legend(ax=ax1, loc=1, frameon=False)

# DISPLAY THE GC% PROFILE ALONG THE SEQUENCE

window_size = 50
windowed_gc_content = [
    gc_content(record.seq[i : i + window_size])
    for i in range(len(record.seq) - window_size)
]
indices = np.arange(len(record.seq) - window_size) + 25

ax2.fill_between(indices[::50], windowed_gc_content[::50], alpha=0.3)
ax2.set_ylabel("GC(%)", fontsize=14)
ax2.set_ylim(bottom=0)

# SAVE THE FIGURE

fig.tight_layout()
fig.savefig("D_display_with_gc_content.svg", bbox_inches="tight")
