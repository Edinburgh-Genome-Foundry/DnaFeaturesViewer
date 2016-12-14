from dna_features_viewer import BiopythonTranslator
from Bio import SeqIO

record = SeqIO.read("example_sequence.gb", "genbank")
graphic_record = BiopythonTranslator().translate_record(record)
ax, _ = graphic_record.plot(figure_width=10)
ax.figure.tight_layout()
ax.figure.savefig("from_genbank.png")
