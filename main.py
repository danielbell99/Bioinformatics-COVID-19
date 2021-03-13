import dataimport
import Nucleotides
import Clustering
import Syntheses
import PairwiseSequencing
import Alignment
import MultipleSequenceAlignment
import test


""" Import Genomes """
# API Parameters
db = 'nucleotide'  # 'protein'
ids = ['844174433']  # append list of IDs you want to fetch from NCBI
genomes = dataimport.read_genomes(db, ids)  # static

# Genome Dictionaries
MERS = genomes[0]
SARSCoV2 = genomes[1]
SARS = genomes[2]
# Create a variable here for a new genome sequence, following the above format e.g.
# NewStrain = genomes[3]


""" Nucleotide Composition """
# AT/GC Content
# Nucleotides.base_content(MERS)  # MERS
# Nucleotides.base_content(SARS)  # SARS
# Nucleotides.base_content(SARSCoV2)  # SARSCoV2

# Polymer Bases Combinations
dimers = Nucleotides.base_combinations(2)  # dinucleotide
# trimers = Nucleotides.base_combinations(3)  # trinucleotide
# tetramers = Nucleotides.base_combinations(4)  # tetranucledotide
# pentamers = Nucleotides.base_combinations(5)  # pentanucledotide
# hexamers = Nucleotides.base_combinations(6)  # hexanucledotide

Nucleotides.composition_comparison(dimers, MERS, SARS, SARSCoV2)
# Nucleotides.composition_comparison(trimers, MERS, SARS, SARSCoV2)
# Nucleotides.composition_comparison(tetramers, MERS, SARS, SARSCoV2)
# Nucleotides.composition_comparison(pentamers, MERS, SARS, SARSCoV2)
# Nucleotides.composition_comparison(hexamers, MERS, SARS, SARSCoV2)


""" Cluster Analysis """
# Clustering.read_normalised_frequencies(dimers, MERS, SARS, SARSCoV2)
# Clustering.read_normalised_frequencies(trimers, MERS, SARS, SARSCoV2)
# Clustering.read_normalised_frequencies(tetramers, MERS, SARS, SARSCoV2)
# Clustering.read_normalised_frequencies(pentamers, MERS, SARS, SARSCoV2)
# Clustering.read_normalised_frequencies(hexamers, MERS, SARS, SARSCoV2)


""" Protein Syntheses """
# Syntheses.protein(MERS)
# Syntheses.protein(SARS)
# Syntheses.protein(SARSCoV2)


""" Pairwise Sequencing """
# PairwiseSequencing.run("DNA", MERS, SARS)
# PairwiseSequencing.run("DNA", MERS, SARSCoV2)
# PairwiseSequencing.run("DNA", SARS, SARSCoV2)

# PairwiseSequencing.run("Protein", MERS, SARS)
# PairwiseSequencing.run("Protein", MERS, SARSCoV2)
# PairwiseSequencing.run("Protein", SARS, SARSCoV2)


""" Alignment Files """
# All sequences considered in Alignment files (.aln)
# Alignment.create_file("DNA")
# Alignment.create_file("Protein")

""" Multiple Sequence Alignment """
""" Bokeh does not support outputs w/ PyCharm (Google for evidence)
Test run Jupyter Notebook version
#MultipleSequenceAlignment.run("DNA")
#MultipleSequenceAlignment.run("Protein")
"""

""" Testing """
# test.run()


print("END OF MAIN")
