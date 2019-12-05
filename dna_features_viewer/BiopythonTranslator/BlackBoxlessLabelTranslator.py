from .BiopythonTranslator import BiopythonTranslator

class BlackBoxlessLabelTranslator(BiopythonTranslator):
    """Translates Biopython records into GraphicRecords where annotations
    appear black on a white background with no box. Which can be cleaner.""" 
    def compute_feature_box_linewidth(self, feature):
        return 0

    def compute_feature_box_color(self, feature):
        return "white"