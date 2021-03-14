import os
from math import log
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import StandardFunctions as sf


def save(bio_type, output_name, sequencing):
    """Stores matches and their scores w/ appropriate file naming convention.

    :param str bio_type: "dna" or "protein"
    :param str output_name: segment of filename that lists virus names (e.g. "_MERS-MT387202_SARS-CoV-JQ316196")
    :param list sequencing: stores all possible matches and scores them
    """
    # Save Pairwise Sequencing
    path_filename = os.path.join('data\\Pairwise Sequencing',
                                 bio_type + output_name + '.txt')  # Defines path and filename
    output = open(path_filename, 'w')
    for a in sequencing:
        output.write(format_alignment(*a))  # Standardised format for output
    output.close()


def gap_function(x, y):
    """Deducts points for gaps when matching up sequencing, using a logarithmic scale.

    :param int x: index to start of gap
    :param int y: length of gap instance
    :return: log() producing negative decimal number
    """
    if y == 0:  # Asks if there's a gap
        return 0
    elif y == 1:  # Gap
        return -1
    # Consecutive gaps
    return round(-((log(y) / 2) + (2 + (y / 4))), 2)  # y = 1 ; return -2.25


def run(bio_type, *names):
    """Performes Pairwise Sequencing task.
    Matches subsequences of similarity that may indicate evolutionary, functional and structural relationships.
    Higher the score, higher the relationship.

    :param bio_type: "dna" or "protein"
    :param ndarray *names: array of str genome names ('*' >=1 genomes)
    """
    if bio_type.lower() == "dna":
        prefix = ""  # none
        directory = 'src/'
    elif bio_type.lower() == "protein":
        prefix = "protein_"
        directory = 'data/Syntheses/'
    else:
        print("biotype: \'" + bio_type + "\' not recongnised. Enter \'DNA\' or \'Protein\'\n")
        return  # Exception Handling

    # Capture all filenames of 'bio_type'
    filenames = sf.capture_filenames(bio_type)
    # filenames = [f.split(".")[0] for f in filenames if bio_type.lower() == "dna"]  # DNA - remove .fasta/.fna and/or .<num> in name

    # Captures two 'bio_type' sequences of interest
    genomes = []  # Genome name, complete protein sequence
    for f in range(len(filenames)):
        for name in names:  # range(len(genomes))
            if set(name).issubset(filenames[f]):
                with open(directory + filenames[f]) as s:
                    # DNA, remove top/ description line
                    seq = ''.join([sf.dna_sequence(sf.output_filename(filenames[f])) if bio_type.lower() == "dna" else s.read()])  # else, protein
                    seq = seq.strip()  # removes any spacing surrounding sequence
                    # Store sequences of interest
                    genome_dict = {'name': name, 'sequence': seq}
                    genomes.append(genome_dict)

    # Concatenates char members of genome as a string object
    seq1 = genomes[0]['sequence']  # Target seq
    seq2 = genomes[1]['sequence']  # Query seq

    # Print Global Sequencing of our pair
    # Match Score - matched sequence chars found; otherwise - Mismatch Score
    # 2 pts, for identical characters
    # -1 pts, for non-identical character
    # -1 pts, when there's a new/ separate gap in the sequence
    # -2.25 pts, if when there's a next gap, right after an existing gap
    sequencing = pairwise2.align.globalmc(seq1[1:100], seq2[1:100], 2, -1, gap_function, gap_function,
                                          penalize_end_gaps=False)
    output_name = sf.output_name([genomes[0]['name'], genomes[1]['name']])  # Appended to 'path_filename', in 'save()'
    save(bio_type, output_name, sequencing)
