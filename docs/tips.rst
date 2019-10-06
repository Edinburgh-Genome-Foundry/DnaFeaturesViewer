Tips
====


Setting the default font
------------------------

The simplest way to change the font for all texts plotted is to change
the Matplotlib default with the following lines at the beginning of the script:

.. code:: python

    import matplotlib
    matplotlib.rcParams['font.family'] = "Walter Turncoat"

.. figure:: _static/images/example_default_font_rc.png
    :align: center

You can also set a default for the label text only by providing a
``default_font_family`` to either a GraphicRecord instance, or (for a global
effect) to the GraphicRecord class itself:

.. code:: python

    graphic_record = GraphicRecord(...)
    graphic_record.default_font_family = 'Walter Turncoat'

    # OR FOR A GLOBAL EFFECT:
    GraphicRecord.default_font_family = 'Walter Turncoat'

.. figure:: _static/images/example_default_font.png
    :align: center