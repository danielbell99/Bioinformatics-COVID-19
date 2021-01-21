def write_file(coronavirus_name, protein):
    """ Stores synthesised protein w/ appropriate file naming convention

    :param str coronavirus_name: name of virus that we've synthesised
    :param protein: actual synthesised data from 'protein()'
    """
    file = open("data\\Syntheses\\" + "protein_" + coronavirus_name, "w")
    file.write(protein)
    file.close()
    print("\n-- Sucessfully Synthesised Protein: " + coronavirus_name + " --\n" + str(protein))

    return


def protein(coronavirus):
    """ DNA -(transcription)> mRNA -(translation)> Protein
    Start codon: AUG (not needed for 'dna_codons' dict)
    Stop codons: UAG, UAA, UGA (denoted by '_')

    :param coronavirus:
    """
    dna = ''.join(coronavirus['sequence'])

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
        'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'TAA': '_', 'TAC': 'Y',
        'TAG': '_', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S',
        'TCT': 'S', 'TGA': '_', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
        'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'
    }
    protein = ""
    if len(dna) % 3 == 2: dna = dna[:-2]  # Ribosome molecular machine only transcribes trinucletides
    if len(dna) % 3 == 1: dna = dna[:-1]  # so ignores any remaining, up to dinucleotide inclusive

    for i in range(0, len(dna), 3):  # increments of 3, i.e. one codon per iter
        codon = dna[i:i + 3]  # a trinucleotide that corresponds to a protein
        protein += dna_codons.get(codon, "")  # appends as str

    write_file(coronavirus['name'], protein)  # Save Synthesised Protein
