import Coronaviridae
import Nucleotides
import Clustering
import Syntheses
import PairwiseSequencing
import Alignment
import MultipleSequenceAlignment
import Testing


""" Import Genomes """
# Genomes
#coronaviridae = Coronaviridae.read_genomes()

# Genome Dictionaries
#MERS = coronaviridae[0]
#SARSCoV2 = coronaviridae[1]
#SARS = coronaviridae[2]
# Create a variable here for a new genome sequence, following the above format e.g.
# NewStrain = coronaviridae[3]


""" Nucleotide Composition """
# AT/GC Content
#Nucleotides.base_content(MERS)  # MERS
#Nucleotides.base_content(SARS)  # SARS
#Nucleotides.base_content(SARSCoV2)  # SARSCoV2

# Polymer Bases Combinations
#dimers = Nucleotides.base_combinations(2)  # dinucleotide
#trimers = Nucleotides.base_combinations(3)  # trinucleotide
#tetramers = Nucleotides.base_combinations(4)  # tetranucledotide

#diComp = Nucleotides.composition_comparison(dimers, MERS, SARS, SARSCoV2)
#triComp = Nucleotides.composition_comparison(trimers, MERS, SARS, SARSCoV2)
#tetraComp = Nucleotides.composition_comparison(tetramers, MERS, SARS, SARSCoV2)


""" Cluster Analysis """
#Clustering.read_normalised_frequencies(trimers, MERS['name'], SARS['name'], SARSCoV2['name'])


""" Protein Syntheses """
#Syntheses.protein(MERS)
#Syntheses.protein(SARS)
#Syntheses.protein(SARSCoV2)


""" Pairwise Sequencing """
#PairwiseSequencing.run("DNA", MERS, SARS)
#PairwiseSequencing.run("DNA", MERS, SARSCoV2)
#PairwiseSequencing.run("DNA", SARS, SARSCoV2)

#PairwiseSequencing.run("Protein", MERS, SARS)
#PairwiseSequencing.run("Protein", MERS, SARSCoV2)
#PairwiseSequencing.run("Protein", SARS, SARSCoV2)


""" Alignment Files """
# All sequences considered in Alignment files (.aln)
#Alignment.create_file("DNA")
#Alignment.create_file("Protein")

""" Multiple Sequence Alignment """
""" Bokeh does not support outputs w/ PyCharm (Google for evidence)
Test run Jupyter Notebook version
#MultipleSequenceAlignment.run("DNA")
#MultipleSequenceAlignment.run("Protein")
"""


""" Testing """
Testing.run()


print("END OF MAIN")