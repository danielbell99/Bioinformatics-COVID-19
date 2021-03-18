import os
from math import log
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.pairwise2 import format_alignment
from joblib.numpy_pickle_utils import xrange
import StandardFunctions as sf


def save(bio_type, output_name, alignments, seq_id, gapless_seq_id):
    """Stores matches and their scores w/ appropriate file naming convention.

    :param str bio_type: "dna" or "protein"
    :param str output_name: segment of filename that lists virus names (e.g. "_MERS-MT387202_SARS-CoV-JQ316196")
    :param Alignment alignments: stores all alignments btwn. a pair of sequences with their match score
    :param float seq_id: match % w/ gaps
    :param float gapless_seq_id: match % w/out gaps
    """
    # File
    path_filename = os.path.join('data/Pairwise Sequencing',
                                 bio_type + output_name + '.txt')  # Defines path and filename
    file = open(path_filename, 'w')

    # Contents
    file.write(format_alignment(*alignments[0]))  # *alignments[0] - best alignment
    file.write('\n' + "Sequence Identity: " + str(round(seq_id, 2)) + "%")
    file.write('\n' + "Gapless Identity: " + str(round(gapless_seq_id, 2)) + "%")

    # Storing all alignments
    """
    for a in sequencing:
        output.write(format_alignment(*a))  # Standardised format for output
    """
    file.close()


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


def sequence_identity(align1, align2):
    """Calculates matches btwn a pair of aligned sequences, as %

    :param str align1: seqA, optimally aligned with align2
    :param str align1: seqB, optimally aligned with align1
    :return: float seq_id: match % w/ gaps
    :return: float gapless_seq_id: match % w/out gaps
    """
    # Sequence ID
    length = len(align1)
    matches = [align1[i] == align2[i] for i in xrange(length)]  # True/ False
    seq_id = (100 * sum(matches)) / length

    print("matches", matches)
    print("seq_id", str(round(seq_id, 2)) + "%")

    # Gap ID
    gapless_length = sum([1 for i in xrange(length) if (align1[i] != '-' and align2[i] != '-')])
    gapless_seq_id = (100 * sum(matches)) / gapless_length  # raw %

    print("gapless_length", gapless_length)
    print("gapless_seq_id", str(round(gapless_seq_id, 2)) + "%")

    return seq_id, gapless_seq_id


def run(bio_type, *names):
    """Performes Pairwise Sequencing task.
    Matches subsequences of similarity that may indicate evolutionary, functional and structural relationships.
    Higher the score, higher the relationship.

    :param bio_type: "dna" or "protein"
    :param ndarray *names: array of str genome names ('*' >=1 genomes)
    """
    bio_type = bio_type.lower()

    if bio_type == "dna":
        prefix = ""  # none
        directory = 'src/'
    elif bio_type == "protein":
        prefix = "protein_"
        directory = 'data/Syntheses/'
    else:
        print("biotype: \'" + bio_type + "\' not recongnised. Enter \'DNA\' or \'Protein\'\n")
        return  # Exception Handling

    # Capture all filenames of 'bio_type'
    filenames = sf.capture_filenames(bio_type)
    # filenames = [f.split(".")[0] for f in filenames if bio_type.lower() == "dna"]  # DNA - remove .fasta/.fna and/or .<num> in name

    # Captures 2 'bio_type' sequences of interest
    genomes = []  # Genome name, complete protein sequence
    for f in range(len(filenames)):
        for name in names:  # range(len(genomes))
            if set(name).issubset(filenames[f]):
                with open(directory + filenames[f]) as s:
                    # DNA, remove top/ description line
                    seq = ''.join([sf.dna_sequence(
                        sf.output_filename(filenames[f])) if bio_type.lower() == "dna" else s.read()])  # else, protein
                    seq = seq.strip()  # removes any spacing surrounding sequence
                    # Store sequences of interest
                    genome_dict = {'name': name, 'sequence': seq}
                    genomes.append(genome_dict)

    # Concatenates char members of genome as str objects
    seqA = genomes[0]['sequence']  # Target seq
    seqB = genomes[1]['sequence']  # Query seq

    # Print Global Sequencing of our pair
    # Match Score - matched sequence chars found; otherwise - Mismatch Score
    # 2 pts, for identical characters
    # -1 pts, for non-identical character
    # -1 pts, when there's a new/ separate gap in the sequence
    # -2.25 pts, if when there's a next gap, right after an existing gap
    alignments = pairwise2.align.globalmc(seqA[1:100], seqB[1:100], 2, -1, gap_function, gap_function,
                                          penalize_end_gaps=False)

    dna_submat = {('A', 'A'): 2, ('A', 'C'): -1, ('A', 'G'): 1, ('A', 'T'): -1,
                  ('C', 'A'): -1, ('C', 'C'): 2, ('C', 'G'): -1, ('C', 'T'): 1,
                  ('G', 'A'): -1, ('G', 'C'): -1, ('G', 'G'): 2, ('G', 'T'): -1,
                  ('T', 'A'): -1, ('T', 'C'): -1, ('T', 'G'): -1, ('T', 'T'): 2}

    #alignments = pairwise2.align.localds(seqA[1:100], seqB[1:100], dna_submat, -2, -1, penalize_end_gaps=False)
    """
    matrix = substitution_matrices.load("BLOSUM62")
    #print(matrix)
    """

    # Model Metrics
    best_aln = alignments[0]  # , one_alignment_only=True)
    align1, align2, score, begin, end = best_aln
    print("best_aln", best_aln)
    print("align1", align1)
    print("align2", align2)
    print("score", score)
    print("begin", begin)
    print("end", end)

    # Sequence Identity
    seq_id, gapless_seq_id = sequence_identity(align1, align2)

    output_name = sf.output_name([genomes[0]['name'], genomes[1]['name']])  # so as to only pass one param
    save(bio_type, output_name, alignments, seq_id, gapless_seq_id)
