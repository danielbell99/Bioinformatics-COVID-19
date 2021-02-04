from os import listdir
from os.path import isfile, join
import StandardFunctions as sf


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
            seq = sf.remove_firstline(directory + s)
            seq = seq.strip()  # removes any spacing surrounding sequence
            print("SEQUENCE: " + str(len(seq)) + "\n" + str(seq))
        else:  # i.e. - bio_type.upper() == "PROTEIN":
            with open(directory + s, 'r') as f:
                seq = f.read()
        seq_lengths.append(len(seq))  # length of each sequence
    print(seq_lengths)

    return seq_lengths


def establish_alignment_file(sequence_filenames, seq_lengths, bio_type, directory):
    """Ensures sequences are of equal length, for MSA.
    Done by appending '_'s to short sequences, no. which calculated by the longest sequence.

    :param list sequence_filenames: list of str; names of files to load in
    :param list seq_lengths: list of int; order is synonymous w/ 'sequence_filenames'
    :param str bio_type: "dna" or "protein"
    :param str directory: folder directory (used for on loading in dna or protein sequences)
    """
    file = open('data/Alignments/' + bio_type.lower() + '.aln', 'w')

    n = -1  # sequence number (tb 0-indexed)
    for s in sequence_filenames:
        n += 1
        with open(directory + s, 'r') as sequence:
            name = sf.output_name_filename(s)  # Coronavirus name

            # Content
            if bio_type.upper() == "DNA":  # ignore top/ description line
                content = sf.remove_firstline(directory + s)
                content = content.strip()  # sometimes there is a space at the end of a DNA sequence

                # print("MAX VALUE IN SEQUENCE: " + str(max(seq_lengths)))
                # print("THIS VALUE IN SEQUENCE: " + str(seq_lengths[n]))
                # print("difference: " + str((max(seq_lengths) - seq_lengths[n])))
                blanks_list = ['_'] * ((max(seq_lengths) - seq_lengths[n]) - (2 + (2 * n)))
                blanks = ''.join(blanks_list)

                content = content + blanks
                content = str.join('', content.splitlines())

            else:  # i.e. - bio_type.upper() == "PROTEIN":
                content = sequence.read()
                blanks_list = ['_'] * (max(seq_lengths) - seq_lengths[n])
                blanks = ''.join(blanks_list)
                content += blanks

            file.write(">" + name + "\n" + content + "\n")  # save
    file.close()

    aln = open('data/Alignments/' + bio_type.lower() + '.aln', 'r')
    print(aln.read())

    return


def create_file(bio_type):
    """Wrap-around function conditions 'bio_type' to create an .aln file of 'bio_type' sequences.
    Invokes: 'sequence_lengths()' & 'establish_alignment_file()'.
    Establishes: 'directory' & 'sequence_filenames', based on 'bio_type' sequences.

    :param str bio_type: "dna" or "protein"
    """
    # Capture either all DNA or all Protein file names
    # Determine which files to load in (case-insensitive)
    sequence_filenames = []
    if bio_type.upper() == "DNA":
        directory = 'src/'
        [sequence_filenames.append(s) for s in listdir(directory) if isfile(join(directory, s))]
        print("DNA of: " + str(sequence_filenames))
    elif bio_type.upper() == "PROTEIN":
        directory = 'data/Syntheses/'
        [sequence_filenames.append(s) for s in listdir(directory) if isfile(join(directory, s))]
        print("Synthesised Proteins of: " + str(sequence_filenames))
    else:
        print("Warning in Alignment.py: biotype \"" + bio_type + "\" not recongnised. \nEnter \"DNA\" or \"Protein\"")
        return  # Exception Handling

    # All sequences must be of the same length
    seq_lengths = sequence_lengths(sequence_filenames, bio_type, directory)
    establish_alignment_file(sequence_filenames, seq_lengths, bio_type, directory)

    return
