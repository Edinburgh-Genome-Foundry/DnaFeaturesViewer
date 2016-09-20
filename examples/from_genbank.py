from dna_features_viewer import GraphicRecord
from Bio import SeqIO
with open("example_sequence.gb", "r") as f:
    record = SeqIO.read(f, "genbank")
graphic_record = GraphicRecord.from_biopython_record(
    record,
    fun_label= lambda feature: feature.qualifiers["label"][0]
)
ax, _ = graphic_record.plot(fig_width=10)
ax.figure.tight_layout()
ax.figure.savefig("from_genbank.png")
