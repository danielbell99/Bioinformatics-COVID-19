# DNA -(transcription)> mRNA -(translation)> Protein
# DNA -> mRNA -> tRNA

def writeFile(coronavirus_name, synthesis, data):
    file = open("data\\Syntheses\\" + synthesis + "_" + coronavirus_name, "w+")
    file.write(data)
    file.close()

    return


def protein(coronavirus):
    # DNA to mRNA
    # Start codon: AUG
    # Stop codons: UAG, UAA, UGA (denoted by '!')
    dna = ''.join(coronavirus['sequence'])

    # dna_codons = collections.defaultdict(lambda: 'Key Not found')
    dna_codons = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""
    if len(dna) % 3 == 2: dna = dna[:-2]  # molecular machine transcribes trinucletides only
    if len(dna) % 3 == 1: dna = dna[:-1]  # so ignores any remainding, up to dinucleotide

    for i in range(0, len(dna), 3):
        codon = dna[i:i + 3]
        protein += dna_codons.get(codon, "")

    print("\n--" + coronavirus['name'] + " | Synthesised Protein --\n" + str(protein))

    writeFile(coronavirus['name'], "protein", protein)  # Save Synthesised Protein


def tRNA(coronavirus):
    # writeFile(coronavirus['name'], "tRNA", tRNA)
    pass
