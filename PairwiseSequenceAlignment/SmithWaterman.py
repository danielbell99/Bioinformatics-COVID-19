from PairwiseSequenceAlignment import tasks
import StandardFunctions as sf


def SmithWaterman(bio_type, *args, **kwargs):
    """Local sequence alignment algorithm.

    :param bio_type: "dna" or "protein"
    :param tuple *args: pair of genomes' str names
    :param dict **kwargs: user-defined point scheme for scoring, and readability
    """
    seqA, seqB, match_points, mismatch_points, gap_points = tasks.preparation(bio_type, *args, **kwargs)
    points_scheme = [match_points, mismatch_points, gap_points]  # deterministic

    m, n = len(seqA), len(seqB)  # let m, n denote the lengths of two sequences

    # Path Traceback - pointer matrix
    pointer_matrix = tasks.empty_matrix(m + 1, n + 1)  # stores traceback path (same dimensions)
    # Dynamic Programming Table - score matrix
    score_matrix = tasks.empty_matrix(m+1, n+1)  # '+1' - far-left col & top row represent indexes of sequences' chars

    max_score = 0  # Initial maximum score obtainable
    # Contents & Declare Pointers
    for i in range(1, m+1):
        for j in range(1, n+1):
            above_score = score_matrix[i][j-1] + gap_points
            left_score = score_matrix[i-1][j] + gap_points
            above_left_score = score_matrix[i-1][j-1] + tasks.match(seqA[i-1], seqB[j-1], points_scheme)  # Assign points to position
            score_matrix[i][j] = max(above_score, left_score, above_left_score, 0)  # cannot be a negative number

            # Update w/ Optimal Traceback Path
            if score_matrix[i][j] == above_score:
                pointer_matrix[i][j] = 1  # 1 - left
            if score_matrix[i][j] == left_score:
                pointer_matrix[i][j] = 2  # 2 - above
            if score_matrix[i][j] == above_left_score:
                pointer_matrix[i][j] = 3  # 3 - above_left
            if score_matrix[i][j] >= max_score:
                max_i, max_j = i, j
                max_score = score_matrix[max_i][max_j];
            if score_matrix[i][j] == 0:  # 0 - End of path
                pointer_matrix[i][j] = 0

    i, j = max_i, max_j  # indices of path starting point
    align1, align2 = '', ''  # initial sequences
    # Traceback - follows pointers
    while pointer_matrix[i][j] != 0:  # 0 - End of path
        if pointer_matrix[i][j] == 1:  # 1 - left
            i -= 1
            align1 += seqA[i]
            align2 += '-'
        elif pointer_matrix[i][j] == 2:  # 2 - above
            j -= 1
            align1 += '-'
            align2 += seqB[j]
        elif pointer_matrix[i][j] == 3:  # 3 - above_left
            i -= 1
            j -= 1
            align1 += seqA[i]
            align2 += seqB[j]


    # Score & Alignment
    align1, signs, align2, str_score = tasks.score_alignment(align1, align2, points_scheme)
    alignments = [align1, signs, align2, str_score]

    # Sequence Identity
    seq_identity = tasks.sequence_identity(align1, align2)  # tuple (seq_id, gapless_seq_id)

    output_name = sf.output_name([args[0][0], args[0][1]])  # so as to only pass one param
    tasks.save("SmithWaterman", bio_type, output_name, alignments, seq_identity)