import Coronaviridae
import Nucleotides
import Clustering
import Syntheses
import PairwiseSequencing
import Alignment

#help(Clustering.read_normalised_frequencies)


# Nitrogenous Bases - of interest
bases = ['A', 'C', 'G', 'T']  # Adenine, Cytosine, Guanine, Thymine

# Polymer Bases Combinations
#dimers = Nucleotides.bases_combinations(bases, 2)  # dinucleotide
trimers = Nucleotides.bases_combinations(bases, 3)  # trinucleotide
#tetramers = Nucleotides.bases_combinations(bases, 4)  # tetranucledotide


# Genomes
coronaviridae = Coronaviridae.read_genomes()

# Genome Dictionaries
MERS = coronaviridae[0]
SARSCoV2 = coronaviridae[1]
SARS = coronaviridae[2]


#Nucleotides.composition(bases, MERS)  # MERS
#Nucleotides.composition(bases, SARS)  # SARS
#Nucleotides.composition(bases, SARSCoV2)  # SARSCoV2

#diComp = Nucleotides.composition_comparison(dimers, MERS, SARS, SARSCoV2)
#triComp = Nucleotides.composition_comparison(trimers, MERS, SARS, SARSCoV2)
#tetraComp = Nucleotides.composition_comparison(tetramers, MERS, SARS, SARSCoV2)


Clustering.read_normalised_frequencies(trimers, MERS['name'], SARS['name'], SARSCoV2['name'])


Syntheses.protein(MERS)
Syntheses.protein(SARS)
Syntheses.protein(SARSCoV2)

#Syntheses.tRNA(SARSCoV2)


#PairwiseSequencing.protein(MERS, SARS)
#PairwiseSequencing.protein(MERS, SARSCoV2)
#PairwiseSequencing.protein(SARS, SARSCoV2)


# All sequences considered in Alignment files
Alignment.create_file("DNA")
Alignment.create_file("Protein")


print("END OF MAIN")