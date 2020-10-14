# https://en.wikipedia.org/wiki/Tandem_repeat
import Coronaviridae
import matplotlib.pyplot as plt
import numpy as np

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


dimers = basesCombination(nucleotide_bases, 2)  # dinucleotide
print("Dimers: " + str(dimers))

trimers = basesCombination(nucleotide_bases, 3)  # trinucleotides
print("Trimers: " + str(trimers))

tetramers = basesCombination(nucleotide_bases, 4)  # tetranucledotides
print("Tetramers: " + str(tetramers))

coronaviridae = Coronaviridae.readInGenomes()


def nucleotideComposition(genome_dict):
    dict = {'A': 'Adenine', 'C': 'Cytosine', 'G': 'Guanine', 'T': 'Thymine'}  # Print name

    print("-- " + genome_dict['coronavirus'] + " --")  # coronavirus
    print(genome_dict['description'])  # description
    sequence = genome_dict['sequence']
    print(sequence)  # sequence

    for n in nucleotide_bases:
        n_count = sequence.count(n)
        n_comp = round(n_count / len(sequence) * 100, 2)
        print(dict[n] + " Composition: " + str(n_comp) + "%")


nucleotideComposition(coronaviridae[0]) # MERS
nucleotideComposition(coronaviridae[1]) # SARS
nucleotideComposition(coronaviridae[2]) # SARSCoV2

# Genome Dictionaries
MERS = coronaviridae[0]
SARS = coronaviridae[1]
SARSCoV2 = coronaviridae[2]

print("************************************")
print(MERS)
print("************************************")
print(SARS)
print("************************************")
print(SARSCoV2)


def compositionComparison(polymers, *genomes):
    # Plots a chart -  Normalised polynucleotide frequences of bases composition, for n genomes of interest

    # polymers - holds basesCombination() output array, x axis
    # *genomes - array holds genome dictionaries, to be plotted

    print(genomes)

    # To display relevant titles
    n = len(polymers[0])  # no. nucleotide units in polymer
    polynucleotide = {1: "Nucleotide", 2: 'Dinucleotide', 3: 'Trinucleotide', 4: 'Tetranucleotide', 5: 'Pentanucleotide'}
    n_coronavirus = "Coronaviridae" if len(genomes) > 1 else "Coronavirus"  # plural or singular

    # Plot onwards
    plt.figure(figsize=(50, 25))

    # X axis
    plt.xticks(rotation=60)
    plt.xlabel(polynucleotide[n], fontsize=30)
    plt.title(polynucleotide[n] + " Composition of " + n_coronavirus)

    # Y axis
    plt.ylabel("Normalised " + polynucleotide[n] + " Frequency", fontsize=30)
    plt.ylim(0, 0.10)
    plt.yticks(np.arange(0, 0.10, 0.01))


    # Content (n genomes)
    for g in genomes:
        nf = normalisedFrequencies(polymers, g['sequence'])
        plt.plot(polymers, nf, linewidth=5, color="red")
        plt.legend(g['coronavirus'])

    # Additions
    plt.grid(True)
    plt.show()
    plt.draw()

    # plt.savefig('polymer_composition.png', format="png")


def normalisedFrequencies(polymers, genome_sequence):
    # Calculates no. appearances a polymer subsequence appears in complete genome sequence, as a percentage
    # returns array of percentages, for each polymer, to be plotted in compositionComparison()

    n = len(polymers[0])  # no. nucleotide units in polymer
    normalisedfreq = []
    for p in polymers:
        count = genome_sequence.count(p)
        nf = (count * n) / len(genome_sequence)
        normalisedfreq.append(nf)

    return normalisedfreq


compositionComparison(trimers, MERS, SARS, SARSCoV2)