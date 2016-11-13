Dna Features Viewer
====================

Dna Features Viewer is a Python library to (wait for it...) visualize DNA
features, e.g. from GenBank or Gff files, using the plotting library Matplotlib:

.. figure:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/by_hand.png
    :align: center

Dna Features Viewer is meant to automatically produce simple and clear plots even
for sequences with lots of overlapping features and long labels.
The plots can be output to many different formats (PNG, JPEG, SVG, PDF), e.g.
for report generation or LIMS interfaces.


License
---------

Dna Features Viewer is an open-source software originally written at the `Edinburgh Genome Foundry
<http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer>`_ under the MIT licence.
Everyone is welcome to contribute !

Installation
--------------

If you have PIP installed, just type in a terminal:

.. code:: python

    (sudo) pip install dna_features_viewer

Dna Features Viewer can be installed by unzipping the source code in one directory and using this command:

.. code:: python

    sudo python setup.py install



Examples of use
----------------


Defining the features by hand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this first example we define features "by hand":

.. code:: python

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

.. figure:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/by_hand.png
    :align: center


If we replace `GraphicRecord` by `CircularGraphicRecord` in the code above we obtain
a circular plot of the construct:

.. figure:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/by_hand_circular.png
    :align: center



Reading the features from a GenBank file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DnaFeaturesViewer plays nice with BioPython. As a result it is super easy to plot the content of a GenBank file:

.. code:: python

    from dna_features_viewer import GraphicRecord
    from Bio import SeqIO
    with open("./plasmid.gb", "r") as f:
        record = SeqIO.read(f, "genbank")
    graphic_record = GraphicRecord.from_biopython_record(record)
    graphic_record.plot(fig_width=10)

.. figure:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/from_genbank.png
    :align: center

Displaying the features along with other plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As it uses Matplotlib, Dna Features Viewer can display the features on top of
other sequences statistics, such as the local GC content:

.. code:: python

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
    _, max_y = graphic_record.plot(ax=ax1m , with_ruler=False)

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

.. figure:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/with_plot.png
    :align: center

Dna Features Viewer is pretty minimal in terms of features but easily extensible since it uses Matplotlib as a backend.

Bonus
~~~~~~

As a bonus, here is what to expect when you feed it with a pathologically annotated Genbank file:

.. figure:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/example_overloaded.png
    :align: center
