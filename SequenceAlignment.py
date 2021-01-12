import os
from math import log
from os import listdir
from os.path import isfile, join
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import StandardFunctions as sf

import pathogenie
#from pybioviz import plotters
#import pybioviz
#from bokeh.io import show


def gap_function(x, y):
    if y == 0:  # No gap
        return 0
    elif y == 1:  # Penalty if gap
        return -2
    return -((2 + (y / 4)) + (log(y) / 2))  # Base case


def protein(*genomes):
    # *genomes - array holds genome dictionaries, to be plotted
    # Singular, Pairwise or Multiple Sequence Alignment

    # Capture all filenames from data/Syntheses (ie. Protein Sequences)
    filenames = []
    [filenames.append(g) for g in listdir('data/Syntheses/') if isfile(join('data/Syntheses/', g))]

    # Captures protein sequences of interest only, for Sequence Alignment
    protein_sequences = []  # Virus name, complete protein sequence
    for i in range(len(genomes)):
        if ("protein_" + genomes[i]['name']) in filenames:
            with open('data/Syntheses/' + "protein_" + genomes[i]['name']) as s:
                sequence = s.read()
                # Store seqeunces of interest
                ps = {'name': genomes[i]['name'], 'sequence': sequence}
                protein_sequences.append(ps)

    # Define n sequences to be aligned
    # concatenates char members of genome as a string object
    seq1 = protein_sequences[0]['sequence']  # Target seq
    seq2 = protein_sequences[1]['sequence']  # Query seq
    
    aln = pathogenie.tools.clustal_alignment(seqs=protein_sequences)
    p = plotters.plot_sequence_alignment(aln)
    show(p)

    # Print Global Alignments of our n sequences
    # Match Score - matched protein chars found; otherwise - Mismatch Score
    # 2, identical characters
    # -1, non-identical character
    # -0.5, when there's a new/ separate gap in the sequence
    # -0.1, if when there's a next gap, right after an existing gap
    alignments = pairwise2.align.globalmc(seq1[1:100], seq2[1:100], 2, -1, gap_function, gap_function,
                                          penalize_end_gaps=False)

    output_name = sf.output_name(genomes)  # Appended to output filename

    # Save Sequence Alignment
    path_filename = os.path.join("data\\Sequence Alignment",
                                 "align" + output_name + ".txt")  # Defines path and filename
    output = open(path_filename, "a")  # Establishes file in directory, for Appending
    for a in alignments:
        output.write(format_alignment(*a))  # Standardised format for output
    output.close()

    return
