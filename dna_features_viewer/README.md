# Code organization

This document walks you trough the Geneblocks code. Please request changes if anything is unclear.

- **GraphicFeature.py** implements a class for defining a *GraphicFeature*, which is an annotation (start, end, strand, label) with graphical properties (color, line width, font family...)

- **GraphicRecord/** implements the *GraphicRecord* class, which can plot a set of *GraphicFeatures* using Matplotlib or Bokeh. To keep file sizes acceptable, many methods are implemented in separate files (*bokeh_plots.py*, *matplotlib_plots.py*) and added to *GraphicRecord* via class mixins.

- **CircularGraphicRecord/** implements the *GraphicRecord* class, which inherits from *GraphicRecord* but draws features circularly using custom Matplotlib patches called "arrow-wedge" (defined in file *ArrowWedge.py*).

- **compute_features_levels.py** implements the algorithm for deciding the levels on which the different features (and annotations) are drawn

- **biotools.py** implements generic biology-related methods (reverse_complement, annotation of Biopython records, etc.)
