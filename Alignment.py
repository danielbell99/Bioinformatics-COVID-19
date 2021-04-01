from os import listdir
from os.path import isfile, join
import StandardFunctions as sf


def dna_alignment(sequence_filenames, seq_lengths, directory):
    """Ensures sequences are of equal length, for MSA.
    Done by appending '-'s to short sequences, no. which calculated by the longest sequence.

    :param list sequence_filenames: list of str; names of files to load in
    :param list seq_lengths: list of int; order is synonymous w/ 'sequence_filenames'
    :param str directory: folder directory (used for on loading in dna or protein sequences)
    """
    file = open('data/Alignments/dna.aln', 'w')

    for idx, val in enumerate(zip(seq_lengths, sequence_filenames)):
        with open(directory + val[1], 'r') as sequence:
            name = sf.output_filename(val[1])  # Genome name
            # Content
            content = sf.dna_sequence(sf.output_filename(val[1]))  # ignore top/ description line
            content = content.strip()  # sometimes there is a space at the end of a DNA sequence
            # print("MAX VALUE IN SEQUENCE: " + str(max(seq_lengths)))
            # print("THIS VALUE IN SEQUENCE: " + str(seq_lengths[i]))
            # print("difference: " + str((max(seq_lengths) - seq_lengths[i])))
            blanks_list = ['-'] * int((max(seq_lengths) - val[0]) - (2 + (2 * idx)))
            blanks = ''.join(blanks_list)
            content = content + blanks
            content = str.join('', content.splitlines())

            file.write(">" + name + "\n" + content + "\n")  # save
    file.close()

    aln = open('data/Alignments/dna.aln', 'r')
    print(aln.read())


def protein_alignment(sequence_filenames, seq_lengths, directory):
    """Same as 'dna_alignment()' but for protein

    :param list sequence_filenames: list of str; names of files to load in
    :param list seq_lengths: list of int; order is synonymous w/ 'sequence_filenames'
    :param str directory: folder directory (used for on loading in dna or protein sequences)
    """
    file = open('data/Alignments/protein.aln', 'w')

    for i, s in zip(seq_lengths, sequence_filenames):
        with open(directory + s, 'r') as sequence:
            name = sf.output_filename(s)  # Genome name
            # Content
            content = sequence.read()
            blanks_list = ['-'] * (max(seq_lengths) - i)
            blanks = ''.join(blanks_list)
            content += blanks

            file.write(">" + name + "\n" + content + "\n")  # save
    file.close()

    aln = open('data/Alignments/protein.aln', 'r')
    print(aln.read())


def sequence_lengths(sequence_filenames, bio_type, directory):
    """Captures lengths of all sequences.

    :param list sequence_filenames: list of str; names of files to load in
    :param str bio_type: "dna" or "protein"
    :param str directory: folder directory (used for on loading in dna or protein sequences)
    :return list seq_lengths: stores ints, lengths of sequences (in order of 'sequence_filenames')
    """
    seq_lengths = []  # in order of 'sequence_filenames'
    for s in sequence_filenames:
        if bio_type.upper() == "DNA":  # ignore top line
            seq = sf.dna_sequence(sf.output_filename(s))
            seq = seq.strip()  # removes any spacing surrounding sequence
            # print("SEQUENCE: " + str(len(seq)) + "\n" + str(seq))
        else:  # i.e. - bio_type.upper() == "PROTEIN":
            with open(directory + s, 'r') as f:
                seq = f.read()
        seq_lengths.append(len(seq))  # length of each sequence
    # print(seq_lengths)

    return seq_lengths


def create(bio_type):
    """Wrap-around function conditions 'bio_type' to create an .aln file of 'bio_type' sequences.
    Invokes: 'sequence_lengths()' & 'establish_alignment_file()'.
    Establishes: 'directory' & 'sequence_filenames', based on 'bio_type' sequences.

    :param str bio_type: "dna" or "protein"
    """
    # Capture either all DNA or all Protein file names
    # Determine which files to load in (case-insensitive)
    if bio_type.lower() != "dna" and bio_type.lower() != "protein":
        print("Warning in Alignment.py: biotype \"" + bio_type + "\" not recongnised. Enter \"DNA\" or \"Protein\"")
        return  # Exception Handling

    directory = sf.directory(bio_type)
    sequence_filenames = [s for s in listdir(directory) if isfile(join(directory, s))]
    print("Alignment of " + bio_type + " Sequences: " + str(sequence_filenames))

    # All sequences must be of the same length
    seq_lengths = sequence_lengths(sequence_filenames, bio_type, directory)
    dna_alignment(sequence_filenames, seq_lengths, directory) if bio_type.lower() == "dna" else protein_alignment(
        sequence_filenames, seq_lengths, directory)
