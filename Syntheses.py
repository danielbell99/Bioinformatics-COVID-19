# DNA -(transcription)> mRNA -(translation)> Protein
# dna to mrna to trna
import Bio
import numpy as np

def writeFile(coronavirus, synthesis, data):
    file = open("data\\Syntheses\\" + synthesis + "_" + coronavirus, "w+")
    file.write(data)
    file.close()
    return

def transcription(coronavirus):
    # DNA to mRNA
    dna = ''.join(coronavirus['sequence'])

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
        'TAC': 'Y', 'TAT': 'Y', 'TAA': 'stop', 'TAG': 'stop',
        'TGC': 'C', 'TGT': 'C', 'TGA': 'stop', 'TGG': 'W',
    }
    protein = ""
    if len(dna) % 3 == 0:
        for i in range(0, len(dna), 3):
            codon = dna[i:i+3]
            protein += dna_codons[codon]
    return protein

    print("\n-- Synthesised Protein --\n" + str(protein))

    writeFile(coronavirus['coronavirus'], "protein", protein)

    return protein


def protein(coronavirus):
    transcription(coronavirus)

    pass


def tRNA(coronavirus):
    mRNA = transcription(coronavirus)

    pass
