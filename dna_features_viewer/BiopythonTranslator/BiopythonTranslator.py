from .BiopythonTranslatorBase import BiopythonTranslatorBase


class BiopythonTranslator(BiopythonTranslatorBase):
    """A translator from SeqRecords to dna_features_viewer GraphicRecord.

    This can be subclassed to create custom "themes" (see the example
    ``custom_biopython_translator.py`` in the docs).

    This class is meant to be customized by subclassing and changing the
    methods (``compute_feature_label``, etc.) and/or the attributes
    (``default_feature_color`` etc).

    Attributes
    ----------

    default_feature_color = "#7245dc"
    graphic_record_parameters
      Dictionary containing keyword arguments that will be passed to the
      (Circular)GraphicRecord constructor.

    ignored_features_types
      A list or tuple of strings indicating all the feature types that should
      always be ignored (i.e. not included in the graphic record) by the
      translator.

    label_fields
      This list of strings provides the order in which the different
      attributes of a Genbank feature will be considered, when automatically
      determining the feature label. For instance if the list is
      ["label", "source", "locus_tag"] and a feature has no label but has a
      "source", the "source" will be displayed in the plots.

    Parameters
    ----------

    features_filters
      List of filters (some_biopython_feature) => True/False.
      Only features passing all the filters are kept.
      This only works if you haven't redefined ``compute_filtered_features``.

    features_properties
      A function (feature)=> properties_dict.
    """

    default_feature_color = "#7245dc"
    ignored_features_types = ()
    label_fields = [
        "label",
        "name",
        "gene",
        "product",
        "source",
        "locus_tag",
        "note",
    ]

    def __init__(self, features_filters=(), features_properties=None):
        self.features_filters = features_filters
        self.features_properties = features_properties

    def compute_feature_color(self, feature):
        """Compute a color for this feature.

        If the feature has a ``color`` qualifier it will be used. Otherwise,
        the classe's ``default_feature_color`` is used.

        To change the behaviour, create a subclass of ``BiopythonTranslator``
        and overwrite this method.
        """
        if "color" in feature.qualifiers:
            color = feature.qualifiers["color"]
            if isinstance(color[0], str):
                return "".join(feature.qualifiers["color"])
            else:
                return color
        else:
            return self.default_feature_color

    def compute_feature_fontdict(self, feature):
        """Compute a font dict for this feature."""
        return None

    def compute_feature_box_linewidth(self, feature):
        """Compute a box_linewidth for this feature."""
        return 0.3

    def compute_feature_box_color(self, feature):
        """Compute a box_color for this feature."""
        return "auto"

    def compute_feature_label_link_color(self, feature):
        """Compute the color of the line linking the label to its feature."""
        return "black"

    def compute_filtered_features(self, features):
        """Return the list of features minus the ignored ones.

        By the method keeps any feature whose type is not in
        ignored_features_types and for which all filter(f) pass.
        """
        return [
            f
            for f in features
            if all([fl(f) for fl in self.features_filters])
            and f.type not in self.ignored_features_types
        ]

    def compute_feature_label(self, feature):
        """Compute the label of the feature."""
        label = feature.type
        for key in self.label_fields:
            if key in feature.qualifiers and len(feature.qualifiers[key]):
                label = feature.qualifiers[key]
                break
        if isinstance(label, list):
            label = "|".join(label)
        return label

    def compute_feature_linewidth(self, feature):
        """Compute the edge width of the feature's arrow/rectangle."""
        return 1.0

    def compute_feature_legend_text(self, feature):
        return None

    def compute_feature_html(self, feature):
        """Gets the 'label' of the feature."""
        return self.compute_feature_label(feature)
