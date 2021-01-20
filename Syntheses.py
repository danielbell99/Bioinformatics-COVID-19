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
    if len(dna) % 3 == 1: dna = dna[:-1]  # so ignores any remaining, up to dinucleotide

    for i in range(0,len(dna), 3):
        codon = dna[i:i+3] # a trinucleotide that corresponds to a protein
        protein += dna_codons.get(codon, "")

    print("\n--" + coronavirus['name'] + " | Synthesised Protein --\n" + str(protein))

    writeFile(coronavirus['name'], "protein", protein)  # Save Synthesised Protein
