"""Generate a 2-plot figure with full sequence on the left,
detail of a sub-segment on the right."""

from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt

record = GraphicRecord(sequence=250*"ATGC", features=[
    GraphicFeature(start=5, end=20, strand=+1, color="#ffd700",
                   label="Small feature"),
    GraphicFeature(start=20, end=500, strand=+1, color="#ffcccc",
                   label="Gene 1 with a very long name"),
    GraphicFeature(start=400, end=700, strand=-1, color="#cffccc",
                   label="Gene 2"),
    GraphicFeature(start=600, end=900, strand=+1, color="#ccccff",
                   label="Gene 3")
])
zoom_start, zoom_end = 398, 428 # coordinates of the "detail"
cropped_record = record.crop((zoom_start, zoom_end))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4))

# PLOT THE WHOLE SEQUENCE

ax1.set_title("Whole sequence", loc='left', weight='bold')
record.plot(ax=ax1)
ax1.fill_between((zoom_start, zoom_end), +1000, -1000, alpha=0.15)

# PLOT THE SEQUENCE DETAILS

cropped_record.plot(ax=ax2)
cropped_record.plot_sequence(ax=ax2)
cropped_record.plot_translation(ax=ax2, location=(408, 423),
                                fontdict={'weight': 'bold'})
ax2.set_title("Sequence detail", loc='left', weight='bold')

fig.savefig('overview_and_detail.png', bbox_inches='tight')
