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
#db = 'nucleotide'  # 'protein'
#ids = ['844174433']  # append list of IDs you want to fetch from NCBI
#genomes = dataimport.read_genomes(db, ids)  # static

# Genome Dictionaries
# MERS = genomes[0]
# SARSCoV2 = genomes[1]
# SARS = genomes[2]
# Create a variable here for a new genome sequence, following the above format e.g.
# NewStrain = genomes[3]


""" Nucleotide Composition """
# AT/GC Content
Nucleotides.base_content("SARS-CoV-JQ316196")
Nucleotides.base_content("MERS-MT387202")
Nucleotides.base_content("SARS-CoV-2-MT873892")


#Nucleotides.base_combinations(2)  # Polymer Bases Combinations - 3, 4, 5, 6

#Nucleotides.composition_comparison(2, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892")
#Nucleotides.composition_comparison(3, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892")
#Nucleotides.composition_comparison(4, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892")
#Nucleotides.composition_comparison(5, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892")
#Nucleotides.composition_comparison(6, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892")


""" Cluster Analysis """
#Clustering.read_normalised_frequencies(2, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892")
#Clustering.read_normalised_frequencies(3, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892")
#Clustering.read_normalised_frequencies(4, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892")
#Clustering.read_normalised_frequencies(5, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892")
#Clustering.read_normalised_frequencies(6, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892")


""" Protein Syntheses """
#Syntheses.protein("SARS-CoV-JQ316196")
#Syntheses.protein("MERS-MT387202")
#Syntheses.protein("SARS-CoV-2-MT873892")


""" Pairwise Sequencing """
#PairwiseSequencing.run("DNA", "SARS-CoV-JQ316196", "MERS-MT387202")
#PairwiseSequencing.run("DNA", "MERS-MT387202", "SARS-CoV-2-MT873892")
#PairwiseSequencing.run("DNA", "SARS-CoV-JQ316196", "SARS-CoV-2-MT873892")

#PairwiseSequencing.run("Protein", "SARS-CoV-JQ316196", "MERS-MT387202")
#PairwiseSequencing.run("Protein", "MERS-MT387202", "SARS-CoV-2-MT873892")
#PairwiseSequencing.run("Protein", "SARS-CoV-JQ316196", "SARS-CoV-2-MT873892")


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
#test.run()


print("END OF MAIN")
