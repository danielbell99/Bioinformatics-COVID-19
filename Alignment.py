from os import listdir
from os.path import isfile, join
import fileinput

import StandardFunctions as sf


def sequence_lengths(sequence_filenames, bio_type, dir):
    # Capture lengths of all sequences
    seq_lengths = []  # in order of 'sequence_filenames'
    for s in sequence_filenames:
        if bio_type.upper() == "DNA":  # ignore top line
            seq = sf.ignore_firstline(dir + s)
            seq = seq.strip()  # sometimes there is a space at the end of a DNA sequence
            print("SEQUENCE: " + str(len(seq)) + "\n" + str(seq))
        else:  # i.e. - bio_type.upper() == "PROTEIN":
            with open(dir + s, "r") as f:
                seq = f.read()
        seq_lengths.append(len(seq))  # length of each sequence
    print(seq_lengths)

    return seq_lengths


def alignment_filename(sequence_filenames, bio_type, dir):
    # Name of Alignment file (.aln)
    str_names = ""
    for s in sequence_filenames:
        # Coronavirus name
        name = sf.output_name_filename(s)
        str_names += "_" + name

    str_alignment_filename = bio_type.lower() + str_names + ".aln"  # Establish filename
    print(str_alignment_filename)

    return str_alignment_filename


def establish_alignment_file(str_alignment_filename, sequence_filenames, seq_lengths, bio_type, dir):
    # Sequences must all be of equal length for MSA
    # Append '_'s to each sequence, number of which calculated by the largest sequence

    file = open("data/Alignments/" + str_alignment_filename, "w")

    n = -1  # sequence number
    for s in sequence_filenames:
        n += 1
        with open(dir + s, "r") as sequence:
            name = sf.output_name_filename(s)  # Coronavirus name

            # Content
            if bio_type.upper() == "DNA":  # ignore top line
                content = sf.ignore_firstline(dir + s)
                content = content.strip()  # sometimes there is a space at the end of a DNA sequence

                #print("MAX VALUE IN SEQUENCE: " + str(max(seq_lengths)))
                #print("THIS VALUE IN SEQUENCE: " + str(seq_lengths[n]))
                #print("difference: " + str((max(seq_lengths) - seq_lengths[n])))
                blanks_list = ['_'] * ((max(seq_lengths) - seq_lengths[n]) - (2 + (2 * n)))
                blanks = ''.join(blanks_list)

                content = content + blanks
                content = str.join("", content.splitlines())

            else:  # i.e. - bio_type.upper() == "PROTEIN":
                content = sequence.read()
                blanks_list = ['_'] * ((max(seq_lengths) - seq_lengths[n]) + 4)  # Additional 4 '_'s, for replacing '_END' postfix
                blanks = ''.join(blanks_list)

                content = content.replace("_END", blanks)

            file.write(">" + name + "\n" + content + "\n")  # save

    file.close()

    aln = open("data/Alignments/" + str_alignment_filename, "r")
    print(aln.read())

    return


def create_file(bio_type):
    """ DESCRIPTION
    Wrap-around function conditions bio_type for creating an .aln file

    :param bio_type: "dna" or "protein"
    :type bio_type: string
    """

    # Capture either all DNA or all Protein filenames
    # Determine no. colours to generate no. characters (case-insensitive)
    sequence_filenames = []
    if bio_type.upper() == "DNA":
        dir = "src/"
        [sequence_filenames.append(s) for s in listdir(dir) if isfile(join(dir, s))]
        print("DNA of: " + str(sequence_filenames))
    elif bio_type.upper() == "PROTEIN":
        dir = "data/Syntheses/"
        [sequence_filenames.append(s) for s in listdir(dir) if isfile(join(dir, s))]
        print("Synthesised Proteins of: " + str(sequence_filenames))
    else:
        print("Warning in Alignment.py: biotype \"" + bio_type + "\" not recongnised. \nEnter \"DNA\" or \"Protein\"")
        return  # Exception Handling

    # All sequences must be of the same length
    seq_lengths = sequence_lengths(sequence_filenames, bio_type, dir)
    # Establish .aln filename under convention
    str_alignment_filename = alignment_filename(sequence_filenames, bio_type, dir)

    establish_alignment_file(str_alignment_filename, sequence_filenames, seq_lengths, bio_type, dir)

    return

