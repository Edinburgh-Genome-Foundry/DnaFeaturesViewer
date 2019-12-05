"""An example with GIF generation at the end. How cool is that!

This example requires the Moviepy library  installed (pip install moviepy).

"""
from Bio import Entrez, SeqIO
from dna_features_viewer import BiopythonTranslator

# DOWNLOAD THE PLASMID's RECORD FROM NCBI
Entrez.email = "zulko@egf.org"
handle = Entrez.efetch(
    db="nucleotide", id=1473096477, rettype="gb", retmode="text"
)
record = SeqIO.read(handle, "genbank")

# CREATE THE GRAPHIC RECORD WITH DNA_FEATURES_VIEWER

color_map = {
    "rep_origin": "yellow",
    "CDS": "orange",
    "regulatory": "red",
    "misc_recomb": "darkblue",
    "misc_feature": "lightblue",
}
translator = BiopythonTranslator(
    features_filters=(lambda f: f.type not in ["gene", "source"],),
    features_properties=lambda f: {"color": color_map.get(f.type, "white")},
)
translator.max_line_length = 15
graphic_record = translator.translate_record(record)
ax, _ = graphic_record.plot(figure_width=8, strand_in_label_threshold=7)
ax.figure.savefig("translator_with_custom_colors.png", bbox_inches="tight")

