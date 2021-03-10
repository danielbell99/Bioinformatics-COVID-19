import numpy as np
import randomcolor
from Bio import AlignIO
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
import panel as pn
pn.extension()  # Displays instances of objects from 'panel' library in Jupyter Notebooks

# Constants - appropriate referencing (e.g. iteration, print, file name, visual)
DNA = ['A', 'C', 'G', 'T', '-']  # Adenine, Cytosine, Guanine, Thymine
PROTEIN = ['I', 'M', 'T', 'N', 'K', 'S', 'R', 'L', 'P', 'H', 'Q', 'R', 'V', 'A', 'D', 'E', 'G', 'S', 'F', 'L', 'Y', 'C',
           'W', '*', '-']


def map_colors(sequences, bio_type):
    """Sets a distinct colour for 'bio_type''s each unique chars in 'sequences'.

    :param list sequences: holds Bio.Seq() for each sequence contents from files (w/out description header in .fasta/.fna)
    :param str bio_type: "dna" or "protein"
    :return list color_set: hexadecimal values as str, for all character in all sequences
    """
    chars = [char for s in list(sequences) for char in s]  # characters in our sequences

    rand_color = randomcolor.RandomColor()

    # Determine no. colours to generate no. characters (case-insensitive)
    if bio_type.lower() == "dna":
        colors = rand_color.generate(count=len(DNA))
        color_map = dict(zip(DNA, colors))
    elif bio_type.lower() == "protein":
        colors = rand_color.generate(count=len(PROTEIN))
        color_map = dict(zip(PROTEIN, colors))
    else:
        return  # run() handles exception

    color_set = [color_map[i] for i in chars]

    return color_set


def visualisation(alignment, bio_type):
    """Multiple Sequence Alignment Viewer.
    Global Viewer - overall of coloured sequences' characters.
    Aligned Sequences Viewer - interactive view for analysis, coloured coded characters & ability to scroll.

    :param AlignIO alignment: Alignment file (.aln) that contains either all DNA or Protein sequences ('data/Alignments/')
    :param str bio_type: "dna" or "protein"
    :return bokeh.models.layouts.Column msa: the Multiple Sequence Alignment interactive visualisations
    """
    sequences = [rec.seq for rec in (alignment)]
    colors = map_colors(sequences, bio_type)
    chars = [char for s in list(sequences) for char in s]
    y_range = [rec.id for rec in alignment]

    num_seqs = len(sequences)
    length = len(sequences[0])  # longest sequence is always first in 'sequences' (Alignment.py)

    # ColumnDataSource
    x = np.arange(1, (length + 1))  # evenly spaced values (so as each character box is the same length and height)
    y = np.arange(0, num_seqs, 1)  # [0, 1, 2, ... len(sequences[0]) + 1]
    x, y = np.meshgrid(x, y)  # returns coordinate matrices of 'x' & 'y' arrays for indexing
    x = x.ravel()  # a 1-D array (faster as no memory is copied)
    y = y.flatten()  # derives a flattened copy of 'yy' array (so as to not modify the returned array)
    recty = (y + 0.5)  # balance color box with char
    # source - crucial part that gathers the sequences' 'chars' w/ their assigned 'colors'
    # into their x/y positions on output
    source = ColumnDataSource(dict(x=x, y=y, recty=recty, text=chars, colors=colors))

    x_range = Range1d(0, (length + 1), bounds='auto')

    # -- Global Viewer --
    glb = figure(plot_width=1000, plot_height=50, x_range=x_range, y_range=y_range,
                 min_border=0, toolbar_location='below')
    rects = Rect(x='x', y='recty', width=1, height=1, fill_color='colors', line_color=None)
    glb.add_glyph(source, rects)  # colours
    glb.yaxis.visible = False

    asv_plot_height = (num_seqs * 25)
    asv_length = 100 if length > 100 else length
    asv_x_range = (0, asv_length)

    # -- Aligned Sequences Viewer --
    # migrate horizontally across sequences
    asv = figure(plot_width=1000, plot_height=asv_plot_height, x_range=asv_x_range, y_range=y_range,
                 min_border=0, toolbar_location='below')
    glyph = Text(x='x', y='y', text_font_size='9pt', text_align='center', text_color='black')
    rects = Rect(x='x', y='recty', width=1, height=1, fill_color='colors', line_color='white', fill_alpha=0.5)
    asv.add_glyph(source, glyph)  # sequence characters
    asv.add_glyph(source, rects)  # colours
    asv.grid.visible = False

    msa = gridplot([[glb], [asv]])

    return msa


def run(bio_type):
    """Procedures for MSA are called here.

    :param str bio_type: "dna" or "protein"
    """
    if bio_type.upper() == "DNA":
        filename = 'dna.aln'
    elif bio_type.upper() == "PROTEIN":
        filename = 'protein.aln'
    else:
        print("Warning in Alignment.py: biotype \"" + bio_type + "\" not recongnised. \nEnter \"DNA\" or \"Protein\"")
        return  # Exception Handling

    alignment = AlignIO.read('data/Alignments/' + filename, 'fasta')
    msa = visualisation(alignment, bio_type)
    pn.pane.Bokeh(msa)

    return
