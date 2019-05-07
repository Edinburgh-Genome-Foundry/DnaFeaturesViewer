from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.PDB.Polypeptide import aa1, aa3

def complement(dna_sequence):
    """Return the complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"TACGGC"``.

    Uses BioPython for speed.
    """
    return str(Seq(dna_sequence).complement())


def reverse_complement(sequence):
    """Return the reverse-complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"GCCGTA"``.

    Uses BioPython for speed.
    """
    return complement(sequence)[::-1]


aa_short_to_long_form_dict = {
    _aa1: _aa3[0] + _aa3[1:].lower()
    for (_aa1, _aa3) in zip(aa1 + '*', aa3 + ['*'])
}
def translate(dna_sequence, long_form=False):
    """Translate the DNA sequence into an amino-acids sequence MLKYQT..."""
    result = str(Seq(dna_sequence).translate())
    if long_form:
        result = [aa_short_to_long_form_dict[aa] for aa in result]
    return result

def extract_translation(sequence, location, long_form=False):
    if len(location) == 3:
        start, end, strand = location
    else:
        start, end = location
        strand = 1
    subsequence = sequence[start: end]
    if strand == -1:
        subsequence = reverse_complement(subsequence)
    translation = translate(subsequence, long_form=long_form)
    if strand == -1:
        translation = translation[::-1]
    return translation


def annotate_biopython_record(seqrecord, location="full",
                              feature_type="misc_feature",
                              margin=0, **qualifiers):
    """Add a feature to a Biopython SeqRecord.

    Parameters
    ----------

    seqrecord
      The biopython seqrecord to be annotated.

    location
      Either (start, end) or (start, end, strand). (strand defaults to +1)

    feature_type
      The type associated with the feature

    margin
      Number of extra bases added on each side of the given location.

    qualifiers
      Dictionnary that will be the Biopython feature's `qualifiers` attribute.
    """
    if location == "full":
        location = (margin, len(seqrecord)-margin)

    strand = location[2] if len(location) == 3 else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type
        )
    )