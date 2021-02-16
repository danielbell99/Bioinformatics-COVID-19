from os import listdir
from os.path import isfile, join

""" Abstract functions used across Application """


def capture_filenames(bio_type):
    """Captures all file names in directory, based on 'bio_type'.

    :param str bio_type: "dna" or "protein"
    :return: list filenames: all file names of 'bio_type'
    """
    filenames = []
    if bio_type.lower() == "dna":
        [filenames.append(d) for d in listdir('src/') if isfile(join('src/', d))]
    elif bio_type.lower() == "protein":
        [filenames.append(p) for p in listdir('data/Syntheses/') if isfile(join('data/Syntheses/', p))]
    else:
        return  # invoker handles exception

    return filenames


def output_name(genomes):
    """Concatenates coronavirus names for an output filename.

    :param dict genomes: contains many entire genome data; 'name' value used
    :return str output: names of genomes of interest
    """
    output = ""
    for i in range(len(genomes)):
        output += '_' + genomes[i]['name']  # separated by '_'
    return output


def output_name_filename(sequence_filename):
    """Concatenates coronavirus names for Alignment file (.aln) setup, for DNA or Protein.
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


def remove_firstline(directory_filename):
    """e.g. DNA datasets have a description line on top, followed by a line break.

    :param str directory_filename: absolute file path
    :return str content: holds the sequence data
    """
    file = open(directory_filename, 'r')
    content = file.readlines()[1:]
    content = ''.join(content)
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


# Not in use
@staticmethod
def remove_prefix(text, prefix):
    """Removes prefix of filenames in data/Syntheses (eg. "protein_<VIRUS>" -> "<VIRUS>").

    :param str text: subject to partition
    :param str prefix: e.g. "nf_"
    :return text: the remainder of str var
    """
    if text.startswith(prefix):
        return text[len(prefix):]
    return text
