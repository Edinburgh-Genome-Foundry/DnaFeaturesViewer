import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator
from Bio import SeqIO
import numpy as np

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 5), sharex=True)

# Parse the genbank file, plot annotations
record = SeqIO.read("example_sequence.gb", "genbank")
translator = BiopythonTranslator()
graphic_record = translator.translate_record(record)
ax, levels = graphic_record.plot()
graphic_record.plot(ax=ax1, with_ruler=False)

# Plot the local GC content
def plot_local_gc_content(record, window_size, ax):
    def gc_content(s):
        return 100.0*len([c for c in s if c in "GC"]) / len(s)
    yy = [gc_content(record.seq[i:i+window_size])
          for i in range(len(record.seq)-window_size)]
    xx = np.arange(len(record.seq)-window_size)+25
    ax.fill_between(xx, yy, alpha=0.3)
    ax.set_ylabel("GC(%)")

plot_local_gc_content(record, window_size=50, ax=ax2)

# Resize the figure to the right height
fig.tight_layout()
fig.savefig("with_plot.png")
