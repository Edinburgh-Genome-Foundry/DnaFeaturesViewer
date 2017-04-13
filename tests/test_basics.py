"""Basic tests to check that the main examples work."""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dna_features_viewer import (BiopythonTranslator, GraphicFeature,
                                 GraphicRecord, CircularGraphicRecord)
from Bio import SeqIO
import numpy as np

example_genbank = os.path.join('tests', "example_sequence.gb")

def test_by_hand(tmpdir):
    """Test building a GraphicRecord "by hand" """
    features = [
        GraphicFeature(start=5, end=20, strand=+1, color="#ffd700",
                       label="Small feature"),
        GraphicFeature(start=20, end=500, strand=+1, color="#ffcccc",
                       label="Gene 1 with a very long name"),
        GraphicFeature(start=400, end=700, strand=-1, color="#cffccc",
                       label="Gene 2"),
        GraphicFeature(start=600, end=900, strand=+1, color="#ccccff",
                       label="Gene 3"),
    ]

    # PLOT AND EXPORT A LINEAR VIEW OF THE CONSTRUCT
    record = GraphicRecord(sequence_length=1000, features=features)
    ax, _ = record.plot(figure_width=5)
    target_file = os.path.join(str(tmpdir), "by_hand.png")
    ax.figure.savefig(target_file)

    # PLOT AND EXPORT A CIRCULAR VIEW OF THE CONSTRUCT
    circular_rec = CircularGraphicRecord(sequence_length=1000,
                                         features=features)
    ax2, _ = circular_rec.plot(figure_width=4)
    ax2.figure.tight_layout()
    target_file = os.path.join(str(tmpdir), "by_hand_circular.png")
    ax2.figure.savefig(target_file, bbox_inches="tight")


def test_from_genbank(tmpdir):
    graphic_record = BiopythonTranslator().translate_record(example_genbank)
    ax, _ = graphic_record.plot(figure_width=10)
    ax.figure.tight_layout()
    target_file = os.path.join(str(tmpdir), "from_genbank.png")
    ax.figure.savefig(target_file)


def test_plot_with_gc_content(tmpdir):

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 4), sharex=True)

    # Parse the genbank file, plot annotations
    record = SeqIO.read(example_genbank, "genbank")
    graphic_record = BiopythonTranslator().translate_record(record)
    ax, levels = graphic_record.plot()
    graphic_record.plot(ax=ax1, with_ruler=False)

    # Plot the local GC content
    def plot_local_gc_content(record, window_size, ax):
        def gc_content(seq):
            return 100.0*len([c for c in seq if c in "GC"]) / len(seq)
        yy = [gc_content(record.seq[i:i+window_size])
              for i in range(len(record.seq)-window_size)]
        xx = np.arange(len(record.seq)-window_size)+25
        ax.fill_between(xx, yy, alpha=0.3)
        ax.set_ylabel("GC(%)")

    plot_local_gc_content(record, window_size=50, ax=ax2)

    # Resize the figure to the right height
    target_file = os.path.join(str(tmpdir), "with_plot.png")
    fig.tight_layout()
    fig.savefig(target_file)
