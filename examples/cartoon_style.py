"""In this example, we draw features XKCD style.

We use Matplotlib's built-in xkcd() function, and a few tweaks:

- We set record.default_box_color to None to prevent annotations
  box drawing.
- We set the record.default_font_family parameter for a nice font
  for annotations
- We set plt.rcParams["font.family"] for a nice font for the
  ruler.
"""
import matplotlib.pyplot as plt
from dna_features_viewer import (GraphicFeature, GraphicRecord,
                                 CircularGraphicRecord)
features = [
    GraphicFeature(start=20, end=500, strand=+1, color="#ffcccc",
                   label="Gene 1 with a name"),
    GraphicFeature(start=400, end=700, strand=-1, color="#cffccc",
                   label="Gene 2"),
    GraphicFeature(start=600, end=900, strand=+1, color="#0000ff",
                   label="Gene 3"),
]

record = GraphicRecord(sequence_length=1000, features=features)
record.default_box_color = None
record.default_font_family = 'Walter Turncoat'

with plt.xkcd():
    plt.rcParams["font.family"] = 'Permanent Marker' # ruler font
    plt.rcParams["xtick.labelsize"] = 'small'
    ax, _ = record.plot(figure_width=5, annotate_inline=False)
    ax.figure.tight_layout()
    ax.figure.savefig("cartoon_style.png", dpi=200)