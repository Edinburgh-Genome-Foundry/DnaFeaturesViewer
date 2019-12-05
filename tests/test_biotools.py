import textwrap
from dna_features_viewer.biotools import (
    extract_graphical_translation,
    reverse_complement,
    find_narrowest_text_wrap,
)


def test_extract_graphical_translation():
    seq = "ATGGACAGAACAATATAA"
    seq1 = "ATGC" + seq + "GTTC"
    seq2 = "ATGC" + reverse_complement(seq) + "GTTC"

    assert extract_graphical_translation(seq1, (4, 22)) == "MDRTI*"
    assert extract_graphical_translation(seq2, (4, 22, -1)) == "MDRTI*"[::-1]


def test_find_narrowest_text_wrap():
    text = "Chloramphenicol resistance marker"
    naive_wrap = textwrap.wrap(text, 30)
    assert (len(naive_wrap) == 2) and len(naive_wrap[0]) == 26
    narrow_wrap = find_narrowest_text_wrap(text, 30)
    lines = narrow_wrap.split("\n")
    assert len(lines) == 2
    assert max(len(l) for l in lines) == 17
