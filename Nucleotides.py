# https://en.wikipedia.org/wiki/Tandem_repeat
import matplotlib.pyplot as plt
import numpy as np
import StandardFunctions

# Constants - appropriate referencing (print, charts, filenames
base_type = {'A': 'Adenine', 'C': 'Cytosine', 'G': 'Guanine', 'T': 'Thymine'}
polynucleotide = {1: "Nucleotide", 2: 'Dinucleotide', 3: 'Trinucleotide', 4: 'Tetranucleotide', 5: 'Pentanucleotide'}
polymer_type = {1: "monomers", 2: "dimers", 3: "trimers", 4: "tetramers", 5: "pentamers"}

def basesCombinations(bases, sub_length):
    n = len(bases)  # no. Bases
    combinations = []  # Appended by basesCombinationRecursive(), externally
    basesCombinationsRecursive(bases, "", n, sub_length, combinations)
    return combinations

def basesCombinationsRecursive(bases, string, n, str_length, combinations):
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
        basesCombinationsRecursive(bases, newString, n, str_length - 1, combinations)

    return combinations


def composition(bases, genome):
    print("\n-- " + genome['name'] + " --")  # name
    print(genome['description'])
    print(genome['sequence'])

    for n in bases:
        n_count = genome['sequence'].count(n)
        n_comp = round(n_count / len(genome['sequence']) * 100, 2)
        print(base_type[n] + " Composition: " + str(n_comp) + "%")

    basesContent(genome['name'], genome['sequence'])

def basesContent(name, sequence):
    # Nitrogenous Bases - AT/GC Ratio:
    # Guanine & Cytosine as % of DNA (always paired together)
    # Adenine & Thymine as % of DNA (always paired together)

    # Data
    gc = ((sequence.count("G") + sequence.count("C")) / len(sequence)) * 100  # GC Content (%)
    at = ((sequence.count("A") + sequence.count("T")) / len(sequence)) * 100  # AT Content (%)
    # ratio = (sequence.count("A") + sequence.count("T")) / (sequence.count("G") + sequence.count("C"))

    # Plot
    plt.title(name + ' Nitrogenous Base Pairs')
    plt.ylabel('Content as Percentage (%)')
    x = ["G-C Content", "A-T Content"]
    y = [gc, at]
    plt.bar(x, y, color='blue')

    # Display & Save
    fig = plt.gcf()
    plt.show()
    plt.draw()
    fig.savefig("data\\Content\\content_" + name + '.png', format="png")


def compositionComparison(polymers, *genomes):
    # Plots a chart -  Normalised polynucleotide frequences of bases composition, for n genomes of interest
    # polymers - holds basesCombination() output array, x axis
    # *genomes - array holds genome dictionaries, to be plotted

    # To display relevant labels
    n = len(polymers[0])  # no. bases in each polymer
    n_coronavirus = "Coronaviridae" if len(genomes) > 1 else "Coronavirus"  # plural or singular
    colormap = ["red", "green", "blue", "yellow", "purple"]

    # Plot
    plt.figure(figsize=(50, 25))
    plt.title(polynucleotide[n] + " Composition of " + n_coronavirus, fontsize=15)

    # X axis
    plt.xlabel(polynucleotide[n], fontsize=10)
    plt.xticks(rotation=45)

    # Y axis
    plt.ylabel("Normalised " + polynucleotide[n] + " Frequency", fontsize=10)
    plt.yticks(np.arange(0, 0.25, 0.01))
    plt.ylim(0, 0.25)

    for g, c in zip(genomes, colormap):
        # Parallel iteration - g (dictionaries in genomes) & c (colours in colourmap)
        nf = normalisedFrequencies(polymers, g)
        plt.plot(polymers, nf, linewidth=2, color=c, label=g['name'])
    plt.legend()

    plt.grid(True)

    # Display & Save
    fig = plt.gcf()
    plt.show()
    plt.draw()

    output_name = StandardFunctions.output_name(genomes)  # Output filename

    fig.savefig("data\\Composition\\" + polynucleotide[n] + output_name + '.png', format="png")

def normalisedFrequencies(polymers, genome):
    # Calculates no. appearances a polymer subsequence appears in complete genome sequence, as a percentage
    # returns array of infinitesimals, for each polymer, to be plotted in compositionComparison()
    # Stores csv - headers = polymers, nf = content
    # filename - "nf_[polymers]_[coronavirus].csv"

    n = len(polymers[0])  # no. nucleotide units in polymer
    sequence = ''.join(genome['sequence'])  # concatenates char members of genome as a string object

    normalisedfreq = []
    for p in polymers:
        count = sequence.count(p)
        nf = (count * n) / len(sequence)
        normalisedfreq.append(nf)

        # Store file - list of polymers & normalised frequency scores

    np.savetxt("data\\Normalised Frequency\\nf_" + polymer_type[n] + "_" + genome['name'] + ".csv", normalisedfreq, delimiter=",")

    return normalisedfreq
