""" Methods used across all Classes
"""

def output_name(genomes):
    # Concatenates coronavirus names for an output filename
    output = ""
    for i in range(len(genomes)):
        output += "_" + genomes[i]['name'] # separated by _
    return output


# Not in use
def removePrefix(text, prefix):
    # removes prefix of filenames in data/Syntheses (eg. "protein_VIRUS" -> "VIRUS")
    if text.startswith(prefix):
        return text[len(prefix):]
    return text