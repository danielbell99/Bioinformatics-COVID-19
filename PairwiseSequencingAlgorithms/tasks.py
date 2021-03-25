import os
from Bio.pairwise2 import format_alignment
from joblib.numpy_pickle_utils import xrange
from math import log
import StandardFunctions as sf

''' Tasks performed by all Pairwise Sequencing Algorithms '''


def save(algorithm_name, bio_type, output_name, alignments, seq_identity):
    """Stores matches and their scores w/ appropriate file naming convention.

    :param str algorithm_name: name of algorithm
    :param str bio_type: "dna" or "protein"
    :param str output_name: segment of filename that lists virus names (e.g. "_MERS-MT387202_SARS-CoV-JQ316196")
    :param Alignment alignments: stores all alignments btwn. a pair of sequences and match score
    :param tuple seq_identity: [0] = 'seq_id' (match % w/ gaps) ; [1] = 'gapless_seq_id' (match % w/out gaps)
    """
    # File
    path_filename = os.path.join('data/Pairwise Sequencing/' + algorithm_name + '/', bio_type.lower() + output_name + '.txt')  # Defines path and filename
    file = open(path_filename, 'w')

    print(type(alignments))
    print(alignments)

    # Contents
    if algorithm_name == "pairwise2.align":
        file.write(format_alignment(*alignments[0]))  # *alignments[0] - best alignment
    else:
        for line in alignments: file.write(str(line) + '\n')  # formatted by 'tasks.score_alignment()'
    file.write('\n' + "Sequence Identity: " + str(round(float(seq_identity[0]), 2)) + '%')
    file.write('\n' + "Gapless Identity: " + str(round(float(seq_identity[1]), 2)) + '%')

    """ Storing all alignments
    for a in sequencing: output.write(format_alignment(*a))  # Standardised format for output
    """
    file.close()


def sequence_identity(align1, align2):
    """Calculates matches btwn a pair of aligned sequences, as %

    :param str align1: optimally aligned 'seqA' for 'seqB''s 'align2'
    :param str align2:  optimally aligned 'seqB' for 'seqA''s 'align1'
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
    print("gapless_seq_id", str(round(gapless_seq_id, 2)) + '%')

    return (seq_id, gapless_seq_id)


def empty_matrix(x, y):
    """Creates empty matrix of determined shape for Dynamic Programming Table or Pointer Matrix.

    :param x: row length (m+1)
    :param y: column length (n+1)
    :return: empty matrix of defined size (0 values)
    """
    matrix = []
    for y_ in range(y):
        matrix.append([])
        for x_ in range(x):
            matrix[-1].append(0)  # 0 values

    return matrix


def match(seqA_c, seqB_c, points_scheme):
    """Determine type of match between aligned sequences' characters.

    :param seqA_c: seqA char
    :param seqB_c: seqB char
    :param list point_scheme: deterministic list of points to score (previously defined)
    :return: float: calculated points
    """
    if seqA_c == seqB_c:
        return points_scheme[0]  # 'match_points'
    elif seqA_c != seqB_c:
        return points_scheme[1]  # 'mismatch_points'
    elif seqA_c == '-' or seqB_c == '-':
        return points_scheme[2]  # 'gap_points'
    else:
        return 0  # Exception Handling


def gap_function(idx, len):  # mandatory pass of 'idx' by 'pairwise2.align's 'x'
    """Deducts points for gaps when matching up sequencing, using a logarithmic scale.

    :param int idx: index to start of gap
    :param int len: gap length
    :return: log(): producing negative decimal number
    """
    if len == 0:  # Asks if there's a gap
        return 0
    elif len == 1:  # Gap
        return -1
    # Consecutive gaps
    return round(-((log(len) / 2) + (2 + (len / 4))), 2)  # len = 1 ; return -2.25


def score_alignment(align1, align2, points_scheme):
    """Return format of Alignment Algorithm for output files.

    :param str align1: optimally aligned 'seqA' for 'seqB''s 'align2'
    :param str align2:  optimally aligned 'seqB' for 'seqA''s 'align1'
    :param list point_scheme: deterministic list of points to score (previously defined)
    :return: str align1, symbols, align2, score: cleansed output vars of alignment algorithm
    """
    score = 0
    signs = ""

    # Reset values
    i, j = 0, 0
    align1 = align1[::-1]  # reverse sequences
    align2 = align2[::-1]

    for i in range(0, len(align1)):
        if align1[i] == align2[i]:  # match
            signs += '|'
            score += match(align1[i], align2[i], points_scheme)
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-':  # mismatch
            score += match(align1[i], align2[i], points_scheme)
            signs += '.'
        elif align1[i] == '-' or align2[i] == '-':  # gap
            signs += ' '
            score += points_scheme[2]  # 'gap_points'

    print(align1 + '\n' + signs + '\n' + align2)
    score = '  Score=' + str(score)
    print(score)

    return align1, signs, align2, score


def preparation(bio_type, *args, **kwargs):
    """Loads in datasets as sequences and establishes points scheme for alignment algorithm.

    :param bio_type: "dna" or "protein"
    :param str name1: name of first genome
    :param str name2: name of second genome
    :param dict **kwargs: user-defined point scheme for scoring, and readability
    :return: str & int: sequences and point scheme, as variables for readability
    """
    bio_type = bio_type.lower()
    if bio_type == "dna":
        seqA = sf.dna_sequence(str(args[0][0]))  # Target sequence
        seqB = sf.dna_sequence(str(args[0][1]))  # Query sequence
    elif bio_type == "protein":
        seqA = sf.protein_sequence(str(args[0][0]))  # Passing names as str
        seqB = sf.protein_sequence(str(args[0][1]))
    else:  # Exception Handling
        return print("biotype: \'" + bio_type + "\' not recongnised. Enter \'DNA\' or \'Protein\'\n")

    seqA = seqA[:100]  # !
    seqB = seqB[:100]

    # Points Scheme - extract or default
    match_points = kwargs.get('match', 2.0)
    mismatch_points = kwargs.get('mismatch', -1.0)
    gap_points = kwargs.get('gap', -5.0)

    return seqA, seqB, match_points, mismatch_points, gap_points