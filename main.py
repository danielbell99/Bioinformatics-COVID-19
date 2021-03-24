import dataimport
import Nucleotides
import Clustering
import Syntheses
import PairwiseSequencing
from PairwiseSequencingAlgorithms import NeedlemanWunsch, SmithWaterman
import Alignment
import MultipleSequenceAlignment
import test


""" Import Genomes """
# API Parameters
db = 'nucleotide'  # 'protein'
ids = ['844174433']  # append list of IDs you want to fetch from NCBI
genomes = dataimport.read_genomes(db, ids)  # static


""" Nucleotide Composition """
# AT/GC Content
#Nucleotides.base_content("SARS-CoV-JQ316196")
#Nucleotides.base_content("MERS-MT387202")
#Nucleotides.base_content("SARS-CoV-2-MT873892")


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
bio_type = "DNA" # "Protein"
#genomes = ["SARS-CoV-JQ316196", "MERS-MT387202"]
#genomes = ["MERS-MT387202", "SARS-CoV-2-MT873892"]
genomes = ["SARS-CoV-JQ316196", "SARS-CoV-2-MT873892"]
points_scheme = {'match': 2.0, 'mismatch': -1.0, 'gap': -5.0}

NeedlemanWunsch.NeedlemanWunsch(bio_type, genomes, **points_scheme)
SmithWaterman.SmithWaterman(bio_type, genomes, **points_scheme)
#PairwiseSequencing.run(bio_type, genomes, **points_scheme)


""" Alignment Files """
# All sequences considered in Alignment files (.aln)
#Alignment.create("DNA")
#Alignment.create("Protein")

""" Multiple Sequence Alignment """
""" Bokeh does not support outputs w/ PyCharm (Google for evidence)
Test run Jupyter Notebook version
#MultipleSequenceAlignment.run("DNA")
#MultipleSequenceAlignment.run("Protein")
"""

""" Testing """
#test.run()


print("END OF MAIN")
