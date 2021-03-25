from PairwiseSequencingAlgorithms import tasks
import StandardFunctions as sf


def NeedlemanWunsch(bio_type, *args, **kwargs):
    """Global sequence alignment algorithm.

    :param bio_type: "dna" or "protein"
    :param tuple *args: pair of genomes' str names
    :param dict **kwargs: user-defined point scheme for scoring, and readability
    """
    seqA, seqB, match_points, mismatch_points, gap_points = tasks.preparation(bio_type, *args, **kwargs)
    points_scheme = [match_points, mismatch_points, gap_points]  # deterministic

    m, n = len(seqA), len(seqB)  # let m, n denote the lengths of two sequences

    # Dynamic Programming Table
    score_matrix = tasks.empty_matrix(m+1, n+1)  # '+1' - far-left col & top row represent indexes of sequences' chars

    # Initial Contents
    for i in range(0, m+1): score_matrix[i][0] = i * gap_points
    for j in range(0, n+1): score_matrix[0][j] = j * gap_points

    for i in range(1, m+1):
        for j in range(1, n+1):
            match_result = score_matrix[i-1][j-1] + tasks.match(seqA[i-1], seqB[j-1], points_scheme)  # Assign points to position
            insertion = score_matrix[i][j-1] + gap_points
            deletion = score_matrix[i-1][j] + gap_points
            score_matrix[i][j] = max(match_result, insertion, deletion)

    align1, align2 = '', ''
    i, j = m, n  # Starting point - bottom-right cell
    # Dynamic Programming Table - Perform Traceback
    while i > 0 and j > 0:  # '>0' - exclude top and far-left indexes
        score = score_matrix[i][j]
        above_score = score_matrix[i][j-1]
        left_score = score_matrix[i-1][j]
        above_left_score = score_matrix[i-1][j-1]

        if score == above_score + gap_points:
            j -= 1
            align1 += '-'
            align2 += seqB[j]
        elif score == left_score + gap_points:
            i -= 1
            align1 += seqB[i]
            align2 += '-'
        elif score == above_left_score + tasks.match(seqA[i-1], seqB[j-1], points_scheme):
            i -= 1
            j -= 1
            align1 += seqA[i]
            align2 += seqB[j]

    # Traceback to top-left index - Gaps
    while i > 0:  # leftward
        i -= 1
        align1 += seqA[i]
        align2 += '-'
    while j > 0:  # topward
        j -= 1
        align1 += '-'
        align2 += seqB[j]


    # Score & Alignment
    align1, signs, align2, str_score = tasks.score_alignment(align1, align2, points_scheme)
    alignments = [align1, signs, align2, str_score]  # naming convention - only one alignment in actuality

    # Sequence Identity
    seq_identity = tasks.sequence_identity(align1, align2)  # tuple (seq_id, gapless_seq_id)

    output_name = sf.output_name([args[0][0], args[0][1]])  # so as to only pass one param
    tasks.save("NeedlemanWunsch", bio_type, output_name, alignments, seq_identity)