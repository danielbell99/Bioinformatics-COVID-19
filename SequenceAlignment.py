from math import log
from os import listdir
from os.path import isfile, join
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# Not in use
def removePrefix(text, prefix):
    # removes prefix of filenames in data/Syntheses (eg. "protein_VIRUS" -> "VIRUS")
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def gap_function(x, y):
    if y == 0:  # No gap
        return 0
    elif y == 1:  # Penalty for gap
        return -2
    return -((2 + (y/4)) + (log(y)/2))  # Base case

def proteinSequenceAlignment(*genomes):
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

    # Print Global Alignments of our n sequences
    # Match Score - matched protein chars found; otherwise - Mismatch Score
    # 2, identical characters
    # -1, non-identical character
    # -0.5, when there's a new/ separate gap in the sequence
    # -0.1, if when there's a next gap, right after an existing gap
    alignments = pairwise2.align.globalmc(seq1, seq2, 2, -1, gap_function, gap_function, penalize_end_gaps = False)

    for a in alignments:
        print(format_alignment(*a))  # standardised format for output

    return

