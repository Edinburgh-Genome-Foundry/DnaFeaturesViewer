from dna_features_viewer import BiopythonTranslator

start, end = 1300, 2700


def feature_properties(f):
    """Fade away all features not overlapping with [start, end]"""
    if f.location.end < start or f.location.start > end:
        return dict(color="white", linecolor="grey", label=None)
    return {}


translator = BiopythonTranslator(features_properties=feature_properties)
graphic_record = translator.translate_record("example_sequence.gb")
ax, _ = graphic_record.plot(figure_width=12, elevate_outline_annotations=True)
ax.fill_between(
    [start, end], -10, 10, facecolor="peachpuff", alpha=0.2, zorder=-1
)

ax.figure.savefig('locally_highlighted_record.png', bbox_inches='tight')
