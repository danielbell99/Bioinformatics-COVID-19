import numpy as np
import randomcolor

from Bio import AlignIO

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

import panel as pn
pn.extension()

dna = ['A', 'C', 'G', 'T', '_']  # Adenine, Cytosine, Guanine, Thymine
protein = ['I', 'M', 'T', 'N', 'K', 'S', 'R', 'L', 'P', 'H', 'Q', 'R', 'V', 'A', 'D', 'E', 'G', 'S', 'F', 'L', 'Y', 'C',
           'W', '_']


def map_colors(seqs, bio_type):
    # Sets distinct color for each 'bio_type''s sequence character
    text = [char for s in list(seqs) for char in s]  # protein characters in our sequences

    rand_color = randomcolor.RandomColor()

    # Determine no. colours to generate no. characters (case-insensitive)
    if bio_type.lower() == "dna":
        colors = rand_color.generate(count=len(dna))
        color_map = dict(zip(dna, colors))
    elif bio_type.lower() == "protein":
        colors = rand_color.generate(count=len(protein))
        color_map = dict(zip(protein, colors))
    else:
        return

    color_set = [color_map[i] for i in text]

    return color_set


def view_alignment(aln, bio_type):
    """Multiple Sequence Alignment Viewer"""

    seqs = [rec.seq for rec in (aln)]
    y_range = [rec.id for rec in aln]
    text = [char for s in list(seqs) for char in s]
    colors = map_colors(seqs, bio_type)

    length = len(seqs[0])  # longest sequence
    num_seqs = len(seqs)

    # ColumnDataSource
    x = np.arange(1, (length + 1))  # [0, 1, 2, ... len(seqs[0]) + 1]
    y = np.arange(0, num_seqs, 1)
    xx, yy = np.meshgrid(x, y)
    gx = xx.ravel()
    gy = yy.flatten()
    recty = (gy + 0.5)

    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = (num_seqs * 25)
    x_range = Range1d(0, (length + 1), bounds='auto')
    if (length > 100):
        asv_length = 100
    else:
        asv_length = length
    asv_x_range = (0, asv_length)
    tools = "xpan, reset, save"

    # Global Sequence View
    glb = figure(title=None, plot_width=1000, plot_height=50, x_range=x_range, y_range=y_range, tools=tools,
                 min_border=0, toolbar_location="below")
    rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors", line_color=None)
    glb.add_glyph(source, rects)
    glb.yaxis.visible = False

    # Aligned Sequences Viewer
    # migrate horizontally across sequences
    asv = figure(plot_width=1000, plot_height=plot_height, x_range=asv_x_range, y_range=y_range, tools=tools,
                 min_border=0, toolbar_location="below")
    glyph = Text(x="x", y="y", text="text", text_align="center", text_color="black", text_font_size="9pt")
    rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors", line_color="white", fill_alpha=0.5)
    asv.add_glyph(source, glyph)  # protein characters
    asv.add_glyph(source, rects)  # colours
    asv.grid.visible = False

    msa = gridplot([[glb], [asv]])
    return msa


def run(bio_type):
    if bio_type.upper() == "DNA":
        file = 'dna_MERS-MT387202.1_SARS-CoV-2-MT873892_SARS-CoV-JQ316196.aln'
    elif bio_type.upper() == "PROTEIN":
        file = 'protein_MERS-MT387202_SARS-CoV-2-MT873892_SARS-CoV-JQ316196.aln'
    else:
        print("Warning in Alignment.py: biotype \"" + bio_type + "\" not recongnised. \nEnter \"DNA\" or \"Protein\"")
        return  # Exception Handling

    aln = AlignIO.read("data/Alignments/" + file, 'fasta')
    msa = view_alignment(aln, bio_type)
    pn.pane.Bokeh(msa)

    return
