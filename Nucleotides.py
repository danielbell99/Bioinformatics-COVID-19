import numpy as np
import matplotlib.pyplot as plt
import StandardFunctions as sf

# Constants - appropriate referencing (e.g. iteration, print, file name, visual)
BASE = ['A', 'C', 'G', 'T']  # Nitrogenous Bases
BASE_NAME = {"A": "Adenine", "C": "Cytosine", "G": "Guanine", "T": "Thymine"}
POLYNUCLEOTIDE = {1: "Nucleotide", 2: "Dinucleotide", 3: "Trinucleotide", 4: "Tetranucleotide", 5: "Pentanucleotide",
                  6: "Hexanuceotide"}
POLYMER = {2: "dimers", 3: "trimers", 4: "tetramers", 5: "pentamers", 6: "hexamers"}  # lower-case for file names
COLOUR_MAP = ["red", "green", "blue", "yellow", "purple"]
Y_AXIS_STOP = {2: 0.25, 3: 0.1, 4: 0.05, 5: 0.025, 6: 0.01}


def base_combinations(polymer_len):
    """Accumulator method for invoking recursive method (memory efficiency).
    Exception Handling: ensures at least dimers are passed.

    :param int polymer_len: length of all polymer instances (e.g. 2 - dimers)
    :return list combinations: all possible polymers, each of a given 'polymer_len'
    """
    if polymer_len <= 1 or not (isinstance(polymer_len, int)):  # checks whole number is above 1
        print("Polynucleotides' length must be 2 to 6")  # context term
        return  # Excepton Handling

    combinations = []  # Appended to by bases_combination_recursive(), externally
    base_combinations_recursive("", polymer_len, combinations)

    # Store file - list of polynucleotides
    np.savetxt('data/Polynucleotides/' + POLYMER[polymer_len], combinations, delimiter=',', fmt='%s')


def base_combinations_recursive(polymer, polymer_len, combinations):
    """Generates all possible polymers possible, each of a 'polymer_len'.
    Each recursion, we append a complete 'polymer' to 'combinations'.
    We pass 'combinations' back to 'base_combinations()' and drop the Stack.
    Stack does not need to be unravelled, once threshold is met, per iteration.

    :param str polymer: an instance of a base combination passed to next recursion and appended to 'combinations'; initially blank ("")
    :param int polymer_len: length of all polymer instances (e.g. 2 - dimers)
    :param list combinations: all possible polymers, each of a given 'polymer_len'
    :return list combinations
    """
    # Base case - run out of character space for each 'polymer'
    if polymer_len == 0:
        combinations.append(polymer)
        # print(polymer)  # Complete polymer instance
        return  # continue back to previous recursion's active for loop iteration

    # Starting with all that begin with 'A' ...
    for i in range(len(BASE)):
        # Next Base appended
        new_polymer = polymer + BASE[i]  # we build upon a new instance
        base_combinations_recursive(new_polymer, (polymer_len - 1), combinations)

    return combinations


def base_content(name):
    """Calculates the composition of each polynucleotide in a given genome as a %.

    :param str name: name of genome
    """
    sequence = sf.dna_sequence(name)

    for n in BASE:
        n_count = sequence.count(n)
        n_comp = round(n_count / len(sequence) * 100, 2)
        print(BASE_NAME[n] + " Composition: " + str(n_comp) + "%")

    bases_content_plot(name, sequence)


def bases_content_plot(name, sequence):
    """Nitrogenous Bases - AT/GC Ratio:
    Guanine & Cytosine as % of DNA (always paired together),
    Adenine & Thymine as % of DNA (always paired together).

    :param str name: name of genome
    :param list sequence: DNA as individual monomer bases
    """
    # Data
    gc = ((sequence.count("G") + sequence.count("C")) / len(sequence)) * 100  # GC Content (%)
    at = ((sequence.count("A") + sequence.count("T")) / len(sequence)) * 100  # AT Content (%)
    # ratio = (sequence.count("A") + sequence.count("T")) / (sequence.count("G") + sequence.count("C"))

    # Plot
    plt.title(name + " Nitrogenous Base Pairs")
    plt.ylabel("Content as Percentage (%)")
    x = ["G-C Content", "A-T Content"]
    y = [gc, at]
    plt.bar(x, y, color='blue')

    # Display & Save
    fig = plt.gcf()
    plt.draw()
    plt.show()
    fig.savefig('data\\Content\\content_' + name + '.png', format='png')


def composition_comparison(polymer_len, *names):
    """Plots a chart - Normalised polynucleotide frequences of bases composition, for n genomes of interest.

    :param int polymer_len: used in 'POLYMER[polymer_len]'
    :param ndarray *names: array of str genome names ('*' >=1 genomes)
    """
    if polymer_len <= 1 or not (isinstance(polymer_len, int)):  # checks whole number is above 1
        print("Polynucleotides' length must be 2 to 6")  # context term
        return  # Excepton Handling

    # holds 'base_combinations()' output array, x axis
    polynucleotides = sf.polynucleotides(polymer_len)

    # To display relevant labels
    polymer_num = len(polynucleotides[0])  # no. bases in each polymer
    n_genomes = "Genomes" if len(names) > 1 else "Genome"  # plural or singular

    # Plot
    x_width = 50 + min(512, (4 ** polymer_num))
    plt.figure(figsize=(x_width, 25))
    plt.title(POLYNUCLEOTIDE[polymer_num] + " Composition of " + n_genomes, fontsize=15)

    # X axis
    polynucleotides = polynucleotides[::polymer_num] if polymer_num > 2 else polynucleotides
    plt.xlabel(POLYNUCLEOTIDE[polymer_num], fontsize=10)
    plt.xticks(rotation=45)

    # Y axis
    y_tick = (0.001 * polymer_num if polymer_num > 2 else 0.01)
    stop = Y_AXIS_STOP[polymer_num]
    plt.ylabel("Normalised " + POLYNUCLEOTIDE[polymer_num] + " Frequency", fontsize=10)
    plt.yticks(np.arange(0, stop, y_tick))
    plt.ylim(0, stop)

    # Parallel iteration
    for name, c in zip(names, COLOUR_MAP):
        sequence = sf.dna_sequence(name)
        nf = normalised_frequencies(polynucleotides, name, sequence)
        plt.plot(polynucleotides, nf, linewidth=2, color=c, label=name)
    plt.legend()

    plt.grid(True)

    # Display & Save
    fig = plt.gcf()
    plt.draw()
    plt.show()

    output_name = sf.output_name(names)  # Appended to output filename
    fig.savefig('data\\Composition\\' + POLYNUCLEOTIDE[polymer_num] + output_name + '.png', format='png')


def normalised_frequencies(polynucleotides, name, sequence):
    """Calculates no. appearances each polymer, as a subsequence, appears in a given genome, as a %.
    Stores .csv - headers = polymers, nf = content (filename "nf_[polymers]_[genome].csv").

    :param list polynucleotides: holds 'base_combinations()' output file, x axis
    :param str name: name of genome
    :param str sequence: DNA of genome
    :return ndarray normalised_freq: infinitesimals, for each polymer, to be plotted in 'composition_comparison()'
    """
    polymer_num = len(polynucleotides[0])  # no. nucleotide units in polymer

    normalised_freq = []
    for p in polynucleotides:
        count = sequence.count(p)
        nf = (count * polymer_num) / len(sequence)
        normalised_freq.append(nf)

    # Store file - normalised frequency scores
    np.savetxt('data/Normalised Frequency/' + POLYMER[polymer_num] + '/nf_' + POLYMER[polymer_num] + '_' + name + '.csv',
               normalised_freq, delimiter=',')

    return normalised_freq
