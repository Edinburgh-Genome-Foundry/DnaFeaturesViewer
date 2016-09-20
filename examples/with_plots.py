import matplotlib.pyplot as plt
from dna_features_viewer import GraphicRecord
from Bio import SeqIO
import numpy as np

fig_width = 10
fig, (ax1, ax2) = plt.subplots(2,1, figsize=(fig_width, 5), sharex=True)

# Parse the genbank file, plot annotations
with open("example_sequence.gb", "r") as f:
    record = SeqIO.read(f, "genbank")
graphic_record = GraphicRecord.from_biopython_record(
    record,
    fun_label= lambda feature: feature.qualifiers["label"][0]
)
_, max_y = graphic_record.plot(ax=ax1, with_ruler=False)

# Plot the local GC content
def plot_local_gc_content(record, window_size, ax):
    gc_content = lambda s: 100.0*len([c for c in s if c in "GC"]) / len(s)
    yy = [gc_content(record.seq[i:i+window_size])
          for i in range(len(record.seq)-window_size)]
    xx = np.arange(len(record.seq)-window_size)+25
    ax.fill_between(xx, yy, alpha=0.3)
    ax.set_ylabel("GC(%)")
plot_local_gc_content(record, window_size=50, ax=ax2)

# Resize the figure to the right height
fig.set_size_inches(fig_width, 2 + 0.4*(max_y+2))
fig.tight_layout()
fig.savefig("with_plot.png")
