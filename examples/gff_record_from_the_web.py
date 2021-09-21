import urllib
from io import StringIO
from dna_features_viewer import BiopythonTranslator

# DEFINE FEATURES ASPECTS


def features_properties(f):
    """Mutations get a red label, other features get a pastel color."""
    label = None
    if f.type == "Mutagenesis":
        label = f.qualifiers["Note"][0]
    color = {
        "Mutagenesis": "firebrick",
        "Active site": "yellow",
        "Beta strand": "lightyellow",
        "Chain": "lightcyan",
        "Helix": "honeydew",
        "Initiator methionine": "white",
        "Metal binding": "lightsteelblue",
        "Turn": "moccasin",
    }.get(f.type, "white")
    return dict(color=color, label=label)


# GET THE RECORD FROM UNIPROT

response = urllib.request.urlopen("https://www.uniprot.org/uniprot/P0A7B8.gff")
record_file = StringIO(response.read().decode())

# TRANSLATE AND PLOT THE RECORD

translator = BiopythonTranslator(features_properties=features_properties)
graphic_record = translator.translate_record(record_file)
ax, _ = graphic_record.plot(
    figure_width=15, max_label_length=100, elevate_outline_annotations=True,
)
ax.set_title("Mutation effects in P0A7B8", fontweight="bold", fontsize=16)
ax.figure.savefig("gff_record_from_the_web.png", bbox_inches="tight")
