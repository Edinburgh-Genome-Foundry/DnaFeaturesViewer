.. raw:: html

    <p align="center">
    <img alt="DNA Features Viewer Logo" title="DNA Features Viewer Logo" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/docs/title.png" width="350">
    </p>

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaFeaturesViewer.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaFeaturesViewer
   :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/DnaFeaturesViewer/badge.svg?branch=master
   :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/DnaFeaturesViewer?branch=master


DNA Features Viewer is a Python library to visualize DNA features, e.g. from GenBank or Gff files:

.. raw:: html

    <p align="center">
    <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/by_hand.png" width="500">
    </p>

Dna Features Viewer is meant to automatically produce simple and clear plots even
for sequences with lots of overlapping features and long labels. It plays well with Matplotlib and Biopython.
The plots can be output to many different formats (PNG, JPEG, SVG, PDF), e.g.
for report generation or LIMS interfaces.


License
---------

Dna Features Viewer is an open-source software originally written at the `Edinburgh Genome Foundry
<http://genomefoundry.org>`_ by `Zulko <https://github.com/Zulko>`_
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

If you intend to use the bokeh features, you need to also install Bokeh and Pandas:

.. code:: python

    (sudo) pip install bokeh pandas


Examples of use
----------------


Basic plots
~~~~~~~~~~~~

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
    record.plot(figure_width=5)

.. raw:: html

    <p align="center">
    <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/by_hand.png" width="500">
    </p>


If we replace `GraphicRecord` by `CircularGraphicRecord` in the code above we obtain
a circular plot of the construct:

.. raw:: html

    <p align="center">
    <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/by_hand_circular.png" width="443">
    </p>

It is also possible to generate interactive (browser-based) plots by using ``plot_with_bokeh`` instead of ``plot``:

.. raw:: html

    <p align="center">
    <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/plot_with_bokeh.png" width="800">
    </p>

Nucleotide sequences, translations, and cropping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DNA features viewer allows to plot nucleotide or amino acid sequences under
the record plot:

.. code:: python

    from dna_features_viewer import GraphicFeature, GraphicRecord

    sequence = "ATGCATGCATGCATGCATGCATGCATGC"
    record = GraphicRecord(sequence, features=[
        GraphicFeature(start=5, end=10, strand=+1, color='#ffcccc'),
        GraphicFeature(start=8, end=15, strand=+1, color='#ccccff')
    ])

    ax, _ = record.plot(figure_width=5)
    record.plot_sequence(ax)
    record.plot_translation(ax, (8, 23), fontdict={'weight': 'bold'})
    ax.figure.savefig('sequence_and_translation.png', bbox_inches='tight')

.. raw:: html

    <p align="center">
    <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/sequence_and_translation.png" width="415">
    </p>

This enables for instance to plot an overview of a sequence along with a detailed detail of a sequence subsegment (`full code <https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/blob/master/examples/overview_and_detail.py>`_)

.. code:: python

    ...
    record.plot(ax=ax1)
    cropped_record = record.crop((zoom_start, zoom_end))
    cropped_record.plot(ax=ax2)
    cropped_record.plot_sequence(ax=ax2)
    cropped_record.plot_translation(ax=ax2, location=(408, 423))

.. raw:: html

    <p align="center">
    <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/overview_and_detail.png" width="900">
    </p>


Reading the features from a GenBank file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DnaFeaturesViewer plays nice with BioPython. As a result it is super easy to plot the content of a Biopython record or directly a GenBank file:

.. code:: python

    from dna_features_viewer import BiopythonTranslator
    graphic_record = BiopythonTranslator().translate_record("my_sequence.gb")
    ax, _ = graphic_record.plot(figure_width=10)

.. raw:: html

    <p align="center">
    <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/from_genbank.png" width="900">
    </p>

The class ``BiopythonTranslator`` determines how the genbank information is transformed into graphical features.
It enables to chose which categories of features to plot, the color of the different features.

Displaying the features along with other plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As it uses Matplotlib, Dna Features Viewer can display the features on top of
other sequences statistics, such as the local GC content:

.. code:: python

    import matplotlib.pyplot as plt
    from dna_features_viewer import BiopythonTranslator
    from Bio import SeqIO
    import numpy as np

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 4), sharex=True)

    # Parse the genbank file, plot annotations
    record = SeqIO.read("example_sequence.gb", "genbank")
    graphic_record = BiopythonTranslator().translate_record(record)
    ax, levels = graphic_record.plot()
    graphic_record.plot(ax=ax1, with_ruler=False)

    # Plot the local GC content
    def plot_local_gc_content(record, window_size, ax):
        gc_content = lambda s: 100.0*len([c for c in s if c in "GC"]) / len(s)
        yy = [gc_content(record.seq[i:i+window_size])
              for i in range(len(record.seq)-window_size)]
        xx = np.arange(len(record.seq)-window_size)+25
        ax.fill_between(xx, yy, alpha=0.3)
        ax.set_ylabel("GC(%)")

    plot_local_gc_content(record, window_size=50, ax=ax2)

    # Resize the figure
    fig.savefig("with_plot.png")


.. raw:: html

    <p align="center">
    <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/with_plot.png" width="800">
    </p>

.. figure::
    :align: center

Custom biopython translators
----------------------------

Dna Features Viewer allows to define "themes" by using custom record translators
instead of the default ``BiopythonTranslator``. Here is an example:

.. code:: python

    from dna_features_viewer import BiopythonTranslator

    class MyCustomTranslator(BiopythonTranslator):
        """Custom translator implementing the following theme:

        - Color terminators in green, CDS in blue, all other features in gold.
        - Do not display features that are restriction sites unless they are BamHI
        - Do not display labels for restriction sites
        - For CDS labels just write "CDS here" instead of the name of the gene.

        """

        def compute_feature_color(self, feature):
            if feature.type == "CDS":
                return "blue"
            elif feature.type == "terminator":
                return "green"
            else:
                return "gold"

        def compute_feature_label(self, feature):
            if feature.type == 'restriction_site':
                return None
            elif feature.type == "CDS":
                return "CDS here"
            else:
                return BiopythonTranslator.compute_feature_label(feature)

        def compute_filtered_features(self, features):
            """Do not display promoters. Just because."""
            return [
                feature for feature in features
                if (feature.type != "restriction_site")
                or ("BamHI" in str(feature.qualifiers.get("label", '')))
            ]


    graphic_record = MyCustomTranslator().translate_record("example_sequence.gb")
    ax, _ = graphic_record.plot(figure_width=10)
    ax.figure.tight_layout()
    ax.figure.savefig("custom_bopython_translator.png")

.. figure:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/master/examples/custom_biopython_translator.png
    :align: center
