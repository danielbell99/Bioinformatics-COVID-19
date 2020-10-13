# https://en.wikipedia.org/wiki/Tandem_repeat
import Coronaviridae

nucleotide_bases = ['A', 'C', 'G', 'T']  # Adenine, Cytosine, Guanine, Thymine


def basesCombination(bases, sub_length):
    n = len(bases)  # no. Nucleotides, of interest
    combinations = []  # Appended by basesCombinationRecursive(), externally
    basesCombinationRecursive(bases, "", n, sub_length, combinations)
    return combinations


def basesCombinationRecursive(bases, string, n, str_length, combinations):
    # Recursively creates all combinations of input Nucleotide Bases, each with a limit of sub_length
    # Base case (as we eventually run out of character space, each string)
    if (str_length == 0):
        combinations.append(string)
        # print(string)  # Complete string
        return

    # Starting with all that begin with 'A' ...
    for i in range(n):
        # Next Base appended
        newString = string + bases[i]  # we build upon a new string
        basesCombinationRecursive(bases, newString, n, str_length - 1, combinations)

    return combinations


dinucleotide = basesCombination(nucleotide_bases, 2)
print("Dinucleotide: " + str(dinucleotide))

trinucleotide = basesCombination(nucleotide_bases, 3)
print("Trinucleotide: " + str(trinucleotide))

tetranucledotide = basesCombination(nucleotide_bases, 4)
print("Tetranucledotide: " + str(tetranucledotide))


coronaviridae = Coronaviridae.readInGenomes()

def nucleotideComposition(genome_dict):
    dict = {'A': 'Adenine', 'C': 'Cytosine', 'G': 'Guanine', 'T': 'Thymine'} # Print name

    print("-- " + genome_dict['coronavirus'] + " --") # coronavirus
    print(genome_dict['description']) # description
    sequence = genome_dict['sequence']
    print(sequence) # sequence

    for n in nucleotide_bases:
        n_count = sequence.count(n)
        n_comp = round(n_count / len(sequence) * 100, 2)
        print(dict[n] + " Composition: " + str(n_comp) + "%")

MERS = nucleotideComposition(coronaviridae[0])
SARS = nucleotideComposition(coronaviridae[0])
SARSCoV2 = nucleotideComposition(coronaviridae[0])

