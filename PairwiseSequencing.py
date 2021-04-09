from PairwiseSequenceAlignment import tasks
from Bio import pairwise2
import StandardFunctions as sf


def run(bio_type, *args, **kwargs):
    """Performs Pairwise Sequencing task.
    Matches subsequences of similarity that may indicate evolutionary, functional and structural relationships.
    Higher the score, higher the relationship.

    :param bio_type: "dna" or "protein"
    :param tuple *args: pair of genomes' str names
    :param dict **kwargs: user-defined point scheme for scoring, and readability
    """
    seqA, seqB, match_points, mismatch_points, gap_points = tasks.preparation(bio_type, *args, **kwargs)

    # Print Global Sequencing of pair
    # Match Score - matched sequence chars found; otherwise - Mismatch Score
    # Defaults:
    # 2.0 pts, for identical characters
    # -1.0 pts, for non-identical character
    # -1.0 pts, when there's a new/ separate gap in the sequence
    # -2.25 pts, if when there's a next gap, right after an existing gap
    alignments = pairwise2.align.globalmc(seqA, seqB, match_points, mismatch_points, tasks.gap_function,
                                          tasks.gap_function, penalize_end_gaps=False)

    # Model Metrics
    align1, align2, score, begin, end = alignments[0]
    print("align1", align1)
    print("align2", align2)
    print("score", score)
    # Indexes of Alignment
    print("begin", begin)
    print("end", end)

    # Sequence Identity
    seq_identity = tasks.sequence_identity(align1, align2)  # tuple (seq_id, gapless_seq_id)

    output_name = sf.output_name([args[0][0], args[0][1]])  # so as to only pass one param
    tasks.save("pairwise2.align", bio_type, output_name, alignments, seq_identity)
