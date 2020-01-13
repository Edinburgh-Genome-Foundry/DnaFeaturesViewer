from Bio import Entrez, SeqIO
from dna_features_viewer import annotate_biopython_record

Entrez.email = "dna_features_viewer@example.com"
handle = Entrez.efetch(
    db="nucleotide", id=1473096477, rettype="gb", retmode="text"
)
record = SeqIO.read(handle, "genbank")
annotate_biopython_record(
    record, location=(40, 1800), feature_type="backbone", label="backbone"
)
record.features = [
    f for f in record.features if f.type not in ["gene", "source"]
]
SeqIO.write(record, "plasmid.gb", "genbank")