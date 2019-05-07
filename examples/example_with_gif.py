"""An example with GIF generation at the end. How cool is that!

This example requires the Moviepy library  installed (pip install moviepy).

"""
from Bio import Entrez, SeqIO
import moviepy.editor as mpe
from moviepy.video.io.bindings import mplfig_to_npimage
import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord

# DOWNLOAD THE PLASMID's RECORD FROM NCBI

handle = Entrez.efetch(db="nucleotide", id=1473096477, rettype="gb",
                       retmode="text")
record = SeqIO.read(handle, "genbank")

# CREATE THE GRAPHIC RECORD WITH DNA_FEATURES_VIEWER

color_map = {
    'rep_origin': 'yellow',
    'CDS': 'orange',
    'regulatory': 'red',
    'misc_recomb': 'darkblue',
    'misc_feature': 'lightblue',
    
}
translator = BiopythonTranslator(
    features_filters=(lambda f: f.type not in ['gene', 'source'],),
    features_properties=lambda f: {'color': color_map.get(f.type, 'white')}
)
translator.max_line_length = 15
graphic_record = translator.translate_record(
    record, record_class=CircularGraphicRecord)
graphic_record.labels_spacing = 30

# ANIMATE INTO A GIF WITH MOVIEPY

duration = 5
def make_frame(t):
    top_nucleotide_index = t * graphic_record.sequence_length / duration
    graphic_record.top_position = top_nucleotide_index
    ax, _ = graphic_record.plot(figure_width=8)
    np_image = mplfig_to_npimage(ax.figure)
    plt.close(ax.figure)
    return np_image

clip = mpe.VideoClip(make_frame, duration=duration)
small_clip = clip.crop(x1=60, x2=-60, y1=60, y2=-100).resize(0.5)
small_clip.write_gif("example_with_gif.gif", fps=15)