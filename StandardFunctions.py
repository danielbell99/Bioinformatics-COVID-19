""" Methods used across all Classes
"""


def output_name(genomes):
    # Concatenates coronavirus names for an output filename
    output = ""
    for i in range(len(genomes)):
        output += "_" + genomes[i]['name']  # separated by _
    return output


def output_name_filename(sequence_filename):
    # Concatenates coronavirus names for Alignment file (.aln) setup
    # One at a time- for loops may have different constructs
    # either DNA or Protein
    if ".fna" in sequence_filename:  # DNA
        name = sequence_filename.replace(".fna", "")
    elif ".fasta" in sequence_filename:
        name = sequence_filename.replace(".fasta", "")
    elif "protein_" in sequence_filename:  # Protein
        name = sequence_filename.replace("protein_", "")

    return name


def ignore_firstline(dir_filename):
    file = open(dir_filename, 'r')
    content = file.readlines()[1:]
    content = ''.join(content)
    file.close()

    return content


# Not in use
@staticmethod
def remove_prefix(text, prefix):
    # removes prefix of filenames in data/Syntheses (eg. "protein_VIRUS" -> "VIRUS")
    if text.startswith(prefix):
        return text[len(prefix):]
    return text
