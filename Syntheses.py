import StandardFunctions as sf

DNA_ALPHABET = ['A', 'C', 'G', 'T']


def write_file(name, protein):
    """Stores synthesised protein w/ appropriate file naming convention.

    :param str name: name of genome
    :param protein: synthesised data from 'protein()'
    """
    file = open('data/Syntheses/' + 'protein_' + name, 'w')
    file.write(protein)
    file.close()
    print("\nSynthesised Protein: " + name + "\n" + str(protein) + "\n")

    return


def alphabet_check(name, sequence, alphabet):
    """Compares all characters/ elements of 'sequence' against its bio_type 'alphabet'.

    :param str name: name of genome
    :param list sequence: holds the sequence data

    :return bool flag: cancel synthesisation
    """
    # Non-DNA Alphabet Character
    # for idx, val in enumerate(sequence)
    if any(c not in alphabet for c in list(sequence)):
        print("Protein Syntheses Halted:")
        print(name + "contains non-DNA character(s)'")
        # User Decision
        print("Do you want to continue synthesisation by skipping over these characters?")
        decided = False
        while not decided:
            decision = input("'Yes' or 'No'")
            if decision.lower() == 'no':
                print("cancelled")
                return 1
            elif decision.lower() == 'yes':
                decided = True
                print("continuing")
            else:
                print("try again...")


def protein(name):
    """DNA -(transcription)> mRNA -(translation)> Protein.
    Start codon: AUG (not needed for 'dna_codons' dict)
    Stop codons: TAA, TAG, TGA (denoted by '*')

    :param str name: name of genome
    """
    sequence = sf.dna_sequence(name)

    # Insufficient
    if len(sequence) < 3:
        print("Protein Syntheses CANCELLED:")
        print("Insufficient Polynucleotides. Synthesisation requires 3 or more.")
        return

    # Non-DNA Characters
    if alphabet_check(name, sequence, DNA_ALPHABET): return

    dna = ''.join(sequence)
    # Ribosome molecular machine only transcribes tri-nucletides
    dna_codons = {
        'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T',
        'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AGA': 'R', 'AGC': 'S',
        'AGG': 'R', 'AGT': 'S', 'ATA': 'I', 'ATC': 'I', 'ATG': 'M',
        'ATT': 'I', 'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CGA': 'R',
        'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L',
        'CTG': 'L', 'CTT': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E',
        'GAT': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GTA': 'V',
        'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'TAA': '*', 'TAC': 'Y',
        'TAG': '*', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S',
        'TCT': 'S', 'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
        'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'
    }

    protein = ""
    for i in range(0, len(dna), 3):  # increments of 3, i.e. one codon per iter
        codon = dna[i:i + 3]  # a tri-nucleotide that corresponds to a protein
        #print("Codon: ", codon)
        #print("Protein: ", dna_codons.get(codon, ''))
        protein += dna_codons.get(codon, '')  # '' - ignore any remaining chars, mono or di-nucleotide

    write_file(name, protein)
