"""This example shows how you can very easily flip a plot horizontally if you
need, using Matplotlib's ax.set_xlim() method."""

from dna_features_viewer import BiopythonTranslator

ax = BiopythonTranslator.quick_class_plot("example_sequence.gb")
x1, x2 = ax.get_xlim()
ax.set_xlim(x2, x1)
ax.figure.tight_layout()
ax.figure.savefig("example_with_inverted_x_axis.png")
