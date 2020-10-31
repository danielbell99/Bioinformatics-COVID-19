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
    seq1 = protein_sequences[0]['sequence']  # concatenates char members of genome as a string object
    seq2 = protein_sequences[1]['sequence']

    # Print Global Alignments of our n sequences
    # Match Score - matched protein chars found; otherwise - Mismatch Score
    alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)

    for a in alignments:
        print(format_alignment(*a))  # standardised format for output

    return

