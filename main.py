import dataimport
import Nucleotides
import Clustering
import Syntheses
import PairwiseSequencing
from PairwiseSequenceAlignment import NeedlemanWunsch, SmithWaterman
import Alignment
import MultipleSequenceAlignment
import test


""" Import Genomes """
# API Parameters
#db = 'nucleotide'  # 'protein'
#ids = ['AP017922']  # append list of IDs you want to fetch from NCBI  # 844174433, AP017922
#genomes = dataimport.read_genomes(db, ids)  # static


""" Nucleotide Composition """
# AT/GC Content
#Nucleotides.base_content("SARS-CoV-JQ316196")
#Nucleotides.base_content("MERS-MT387202")
#Nucleotides.base_content("SARS-CoV-2-MT873892")
#Nucleotides.base_content("StaphylococcusAureus-AP017922")

#Nucleotides.base_combinations(2)  # Polymer Bases Combinations - 3, 4, 5, 6

#Nucleotides.composition_comparison(2, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892", "StaphylococcusAureus-AP017922")
#Nucleotides.composition_comparison(3, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892", "StaphylococcusAureus-AP017922")
#Nucleotides.composition_comparison(4, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892", "StaphylococcusAureus-AP017922")
#Nucleotides.composition_comparison(5, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892", "StaphylococcusAureus-AP017922")
#Nucleotides.composition_comparison(6, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892", "StaphylococcusAureus-AP017922")


""" Cluster Analysis """
#Clustering.run(100, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892")  # 'n_reads'
#Clustering.run(100, "SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892", "StaphylococcusAureus-AP017922")


""" Protein Syntheses """
#Syntheses.protein("SARS-CoV-JQ316196")
#Syntheses.protein("MERS-MT387202")
#Syntheses.protein("SARS-CoV-2-MT873892")
#Syntheses.protein("StaphylococcusAureus-AP017922")


""" Pairwise Sequencing """
#bio_type = "DNA"  # "Protein"
#genomes = ["SARS-CoV-JQ316196", "MERS-MT387202"]
#genomes = ["MERS-MT387202", "SARS-CoV-2-MT873892"]
#genomes = ["SARS-CoV-JQ316196", "SARS-CoV-2-MT873892"]

nw_points_scheme = {'match': 1.0, 'mismatch': -1.0, 'gap': -1.0}  # wiki
sw_points_scheme = {'match': 3.0, 'mismatch': -3.0, 'gap': -2.0}  # wiki
ps_points_scheme = {'match': 2.0, 'mismatch': -1.0, 'gap': -5.0}  # Biopython documentation

'''
SmithWaterman.SmithWaterman("DNA", ["SARS-CoV-JQ316196", "MERS-MT387202"], **sw_points_scheme)
print("done")
SmithWaterman.SmithWaterman("DNA", ["MERS-MT387202", "SARS-CoV-2-MT873892"], **sw_points_scheme)
print("done")
SmithWaterman.SmithWaterman("DNA", ["SARS-CoV-JQ316196", "SARS-CoV-2-MT873892"], **sw_points_scheme)
print("done")

SmithWaterman.SmithWaterman("Protein", ["SARS-CoV-JQ316196", "MERS-MT387202"], **sw_points_scheme)
print("done")
SmithWaterman.SmithWaterman("Protein", ["MERS-MT387202", "SARS-CoV-2-MT873892"], **sw_points_scheme)
print("done")
SmithWaterman.SmithWaterman("Protein", ["SARS-CoV-JQ316196", "SARS-CoV-2-MT873892"], **sw_points_scheme)
print("done")
'''

#NeedlemanWunsch.NeedlemanWunsch(bio_type, genomes, **nw_points_scheme)
#SmithWaterman.SmithWaterman(bio_type, genomes, **sw_points_scheme)
#PairwiseSequencing.run(bio_type, genomes, **ps_points_scheme)


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
