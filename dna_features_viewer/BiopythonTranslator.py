from .GraphicRecord import GraphicRecord
from .GraphicFeature import GraphicFeature
from Bio import SeqIO

class BiopythonTranslator:
    """A translator from SeqRecords to dna_features_viewer GraphicRecord.

    Parameters
    ----------

    features_filters
      List of filters (some_biopython_feature) => True/False.
      Only features passing all the filters are kept.

    features_properties
      A function (feature)=> properties_dict

    """
    default_feature_color = "#7245dc"

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
        return feature.qualifiers.get("color", self.default_feature_color)

    def compute_filtered_features(self, features):
        return [f for f in features
                if all([fl(f) for fl in self.features_filters])]

    @staticmethod
    def compute_feature_label(feature):
        """Gets the 'label' of the feature.

        This method looks for the first non-empty qualifier of the feature
        in this order ``label``, ``source``, ``locus_tag``, ``note``, which
        means that you can provide a label with, for instance
        ``feature.qualifiers['note'] = 'some note'``.

        To change the behaviour, create a subclass of ``BiopythonTranslator``
        and overwrite this method.
        """
        result = feature.type
        for key in ["label", "source", "locus_tag", "note"]:
            if key in feature.qualifiers:
                result = feature.qualifiers[key]
                break
        if isinstance(result, list):
            return "|".join(result)
        else:
            return result

    @staticmethod
    def compute_feature_html(feature):
        """Gets the 'label' of the feature."""
        result = feature.type
        for key in ["note", "locus_tag", "label", "source"]:
            if key in feature.qualifiers:
                result = feature.qualifiers[key]

                break
        if isinstance(result, list):
            return "|".join(result)
        else:
            return result


    def translate_feature(self, feature):
        """Translate a Biopython feature into a Dna Features Viewer feature."""
        properties = dict(label=self.compute_feature_label(feature),
                          color=self.compute_feature_color(feature),
                          html=self.compute_feature_html(feature))
        if self.features_properties is not None:
            other_properties = self.features_properties(feature)
        else:
            other_properties = {}
        properties.update(other_properties)

        return GraphicFeature(start=feature.location.start,
                              end=feature.location.end,
                              strand=feature.location.strand,
                              **properties)

    def translate_record(self, record, grecord_class=None):
        """Create a new GraphicRecord from a BioPython Record object.

        Parameters
        ----------

        record
          A BioPython Record object or the path to a Genbank file.

        grecord_class
          The graphic record class to use, e.g. GraphicRecord (default) or
          CircularGraphicRecord.
        """

        if isinstance(record, str):
            record = SeqIO.read(record, "genbank")
        if grecord_class is None:
            grecord_class = GraphicRecord
        return grecord_class(sequence_length=len(record.seq), features=[
            self.translate_feature(feature)
            for feature in self.compute_filtered_features(record.features)
            if feature.location is not None
        ])
