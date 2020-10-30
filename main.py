import Coronaviridae
import Nucleotides
import Clustering
import Syntheses

# Nitrogenous Bases - of interest
bases = ['A', 'C', 'G', 'T']  # Adenine, Cytosine, Guanine, Thymine

# Polymer Bases Combinations
dimers = Nucleotides.basesCombinations(bases, 2)  # dinucleotide
#trimers = Nucleotides.basesCombinations(bases, 3)  # trinucleotide
#tetramers = Nucleotides.basesCombinations(bases, 4)  # tetranucledotide


# Genomes
coronaviridae = Coronaviridae.readInGenomes()

# Genome Dictionaries
MERS = coronaviridae[0]
SARSCoV2 = coronaviridae[1]
SARS = coronaviridae[2]


#Nucleotides.composition(bases, MERS)  # MERS
#Nucleotides.composition(bases, SARS)  # SARS
#Nucleotides.composition(bases, SARSCoV2)  # SARSCoV2

#diComp = Nucleotides.compositionComparison(dimers, MERS, SARS, SARSCoV2)
#triComp = Nucleotides.compositionComparison(trimers, MERS, SARS, SARSCoV2)
#tetraComp = Nucleotides.compositionComparison(tetramers, MERS, SARS, SARSCoV2)


#Clustering.readnNormalisedFrequencies(dimers, MERS['coronavirus'], SARS['coronavirus'], SARSCoV2['coronavirus'])


Syntheses.protein(MERS)
#Syntheses.protein(SARS)
#Syntheses.protein(SARSCoV2)

#Syntheses.tRNA(SARSCoV2)