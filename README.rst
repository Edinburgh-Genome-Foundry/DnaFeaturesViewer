Dna Features Viewer
====================

Dna Features Viewer is a Python library to (wait for it...) visualize DNA
features, e.g. from GenBank or Gff files, using the plotting library Matplotlib:

Dna Features Viewer is fairly minimal (<200 lines of code) but can display sequences with lots
of overlapping features and long labels, without getting too messy.
It also plays nicely with Biopython.


Installation
--------------

Dna Features Viewer can be installed by unzipping the source code in one directory and using this command: ::

    sudo python setup.py install

PIP install is coming soon !

Contribute
-----------

Dna Features Viewer is an open-source software originally written by Zulko
at the Edinburgh Genome Foundry and released under the MIT licence.
Everyone is welcome to contribute !

Examples of use
----------------

In this first example we define features "by hand":
::
    from dna_features_viewer import GraphicFeature, GraphicRecord
    features=[
        GraphicFeature(start=0, end=20, strand=+1, color="#ffd700",
                       label="Small feature"),
        GraphicFeature(start=20, end=500, strand=+1, color="#ffcccc",
                       label="Gene 1 with a very long name"),
        GraphicFeature(start=400, end=700, strand=-1, color="#cffccc",
                       label="Gene 2"),
        GraphicFeature(start=600, end=900, strand=+1, color="#ccccff",
                       label="Gene 3")
    ]
    record = GraphicRecord(sequence_length=1000, features=features)
    record.plot(fig_width=5)


Here is how you parse a GenBank file using BioPython and display the features
using Dna Features Viewer:
::
    from dna_features_viewer import GraphicRecord
    from Bio import SeqIO
    with open("./plasmid.gb", "r") as f:
        record = SeqIO.read(f, "genbank")
    graphic_record = GraphicRecord.from_biopython_record(record)
    graphic_record.plot(fig_width=10)

As it uses Matplotlib, Dna Features Viewer can display the features on top of
other sequences statistics, such as the local GC content:
::
    import matplotlib.pyplot as plt
    from dna_features_viewer import GraphicRecord
    from Bio import SeqIO
    import numpy as np

    figure_width = 10
    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(figure_width,5), sharex=True)

    # Parse the genbank file, plot annotations
    with open("./plasmid.gb", "r") as f:
        record = SeqIO.read(f, "genbank")
    graphic_record = GraphicRecord.from_biopython_record(record)
    _, max_y = graphic_record.plot(ax=ax1)

    # Plot the local GC content
    def plot_local_gc_content(record, window_size, ax):
        gc_content = lambda s: 1.0*len([c for c in s if c in "GC"]) / len(s)
        yy = [gc_content(record.seq[i:i+window_size])
              for i in range(len(record.seq)-window_size)]
        xx = np.arange(len(record.seq)-window_size)+25
        ax.fill_between(xx, yy, alpha=0.3)
    plot_local_gc_content(record, window_size=50, ax=ax2)

    # Resize the figure
    fig.set_size_inches(figure_width, 2 + 0.4*(max_y+2))
