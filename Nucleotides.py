# https://en.wikipedia.org/wiki/Tandem_repeat

nucleotide_bases = ['A', 'C', 'G', 'T']  # Adenine, Cytosine, Guanine, Thymine


def basesCombination(bases, sub_length):
    n = len(bases)  # no. Nucleotides
    combinations = []  # Appended by basesCombinationRecursive(), externally
    basesCombinationRecursive(bases, "", n, sub_length, combinations)
    return combinations

def basesCombinationRecursive(bases, string, n, str_length, combinations):
    # Recursively creates all combinations of input Nucleotide Bases, each with a limit of sub_length
    # Base case (as we eventually run out of character space, each string)
    if (str_length == 0):
        combinations.append(string)
        print(string)  # Complete string
        return

    # Starting with all that begin with 'A' ...
    for i in range(n):
        # Next Base appended
        newString = string + bases[i]  # we build upon a new string
        basesCombinationRecursive(bases, newString, n, str_length - 1, combinations)

    return combinations


dinucleotide = basesCombination(nucleotide_bases, 2)

trinucleotide = basesCombination(nucleotide_bases, 3)

tetranucledotide = basesCombination(nucleotide_bases, 4)
