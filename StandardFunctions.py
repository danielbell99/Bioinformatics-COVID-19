from os import listdir
from os.path import isfile, join

''' Abstract functions used across Application '''


def capture_filenames(bio_type):
    """Captures all file names in directory, based on 'bio_type'.

    :param str bio_type: "dna" or "protein"
    :return: list filenames: all file names of 'bio_type'
    """
    if bio_type.lower() == "dna":
        filenames = [d for d in listdir('src/') if isfile(join('src/', d))]
    elif bio_type.lower() == "protein":
        filenames = [p for p in listdir('data/Syntheses/') if isfile(join('data/Syntheses/', p))]
    else:
        return  # invoker handles exception

    return filenames


def output_name(names):
    """Concatenates genome names for an output filename.

    :param list names: names of genomes
    :return str output: concatenated names of genomes of interest
    """
    output = ""
    for name in names:
        output += '_' + name  # separated by '_'
    return output


def output_filename(sequence_filename):
    """Concatenates genome names for Alignment file (.aln) setup, for DNA or Protein.
    One at a time - for loops may have different constructs.

    :param str sequence_filename: subject to partition
    :return str name: remainder of str var
    """
    if '.fna' in sequence_filename:  # DNA
        name = sequence_filename.replace('.fna', '')
    elif '.fasta' in sequence_filename:
        name = sequence_filename.replace('.fasta', '')
    elif 'protein_' in sequence_filename:  # Protein
        name = sequence_filename.replace('protein_', '')

    return name


def dna_sequence(filename):
    """e.g. DNA datasets have a description line on top, followed by a line break.

    :param str filename: one filename of interest, w/out extension
    :return str content: holds the sequence data
    """
    try:
        file = open('src/' + filename + '.fasta', 'r')
    except:
        file = open('src/' + filename + '.fna', 'r')
    content = file.readlines()[1:]
    content = ''.join(content)
    content = ''.join(content.splitlines())  # bc splitlines() removes '\n' but returns list
    file.close()

    return content


def protein_sequence(filename):
    """Retrieves Protein sequence as str.

    :param str filename: one filename of interest, w/out extension
    :return str content: holds the sequence data
    """
    try:
        file = open('data/Syntheses/protein_' + filename + '.txt', 'r')
    except:
        file = open('data/Syntheses/protein_' + filename, 'r')
    content = file.readlines()
    content = ''.join(content)
    content = ''.join(content.splitlines())  # bc splitlines() removes '\n' but returns list
    file.close()

    return content


def directory(bio_type):
    """Determines parent directory, based on 'bio_type'.

    :param str bio_type: "dna" or "protein"
    :return: str directory: 'src/' or 'data/Syntheses/'
    """
    if bio_type.lower() == "dna":
        directory = 'src/'
    elif bio_type.lower() == "protein":
        directory = 'data/Syntheses/'
    else:
        return  # invoker handles exception

    return directory


def polynucleotides(polymer_len):
    """Returns list of polynucleotides of interest

    :param int polymer_len: used in 'POLYMER[polymer_len]'
    :return str output: concatenated names of genomes of interest
    """
    POLYMER = {2: "dimers", 3: "trimers", 4: "tetramers", 5: "pentamers", 6: "hexamers"}  # lower-case for file names

    file = open('data/Polynucleotides/' + POLYMER[polymer_len], 'r')
    polynucleotides = file.read().splitlines()
    file.close()

    return polynucleotides


# Not in use
# @staticmethod
def remove_prefix(text, prefix):
    """Removes prefix of filenames in data/Syntheses (eg. "protein_<VIRUS>" -> "<VIRUS>").

    :param str text: subject to partition
    :param str prefix: e.g. "nf_"
    :return text: the remainder of str var
    """
    if text.startswith(prefix):
        return text[len(prefix):]
    return text
