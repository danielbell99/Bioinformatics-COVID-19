import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import time
import seaborn as sns
import matplotlib.pyplot as plt
import StandardFunctions as sf

# Constants - appropriate referencing (e.g. iteration, print, file name, visual)
POLYMER = {1: "monomers", 2: "dimers", 3: "trimers", 4: "tetramers", 5: "pentamers",
           6: "hexamers"}  # lower-case for file names
PATH = 'data/Normalised Frequency/'
COLOUR_MAP = ["red", "green", "blue", "yellow", "purple"]


def read_normalised_frequencies(polymer_len, *names):
    """Imports produced Normalised Frequency files (depending on: polynucleotides & genomes of interest).
    Normalised Frequencies are passed through Dimentionality Reduction algorithms for Cluster Analysis.
    Measuring time elapsed.
    Minimum of 2 Genome names required - '*others' is optional.

    :param int polymer_len: used in 'POLYMER[polymer_len]'
    :param ndarray *names: array of str genome names ('*' >=1 genomes)
    """
    genome_names = [name for name in names]
    polymer_name = POLYMER[polymer_len]  # e.g. "dimers" if 2

    # Capture filenames from data/Normalised Frequency
    # Append each file as a column
    nf_df = pd.DataFrame()  # DataFrame passed to PCA (dimensionality reduction algorithm)
    for file, cv_name in zip(os.listdir(PATH + polymer_name + '/'), genome_names):
        file_dir = PATH + polymer_name + '/' + 'nf_' + polymer_name + '_' + cv_name + '.csv'
        print(file_dir)
        if os.path.isfile(os.path.join(file_dir)):
            column = pd.read_csv(os.path.join(file_dir), header=None)
            nf_df = pd.concat([nf_df, column], axis=1)
    # nf_df = nf_df.T  # transpose - rows = viruses, cols = normalised freqs.
    print(nf_df)

    principal_component_analysis(nf_df, genome_names, polymer_len)
    tSNE(nf_df, genome_names, polymer_len)


def seaborn_scatterplot(model, results, genome_names, polymer_len):
    """Standardised method for using Seaborn's Scatterplot to visualise Unsupervised ML results.

    :param str model: String literal name of model for saving figure as filename
    :param ndarray results: contains entire genome data; 'name' value used (*others is optional)
    :param list polynucleotides: combination of nucleic acids of n length each
    :param list genome_names: names extracted from included genomes['name']
    :param int polymer_name: len(polynucleotides[0]) used to get polymer name from POLYMER dict
    """
    str_names = '_'.join(genome_names)
    polynucleotides = sf.polynucleotides(polymer_len)
    polymer_name = POLYMER[polymer_len]  # e.g. "dimers" if 2

    # df = pd.DataFrame(data=results)
    # df = pd.concat([df, polynucleotides], axis=1)
    # results = df.to_numpy()
    # print("HELLO")
    # print(results)

    sns.scatterplot(x=results[:, 0], y=results[:, 1], alpha=1, s=100).plot()
    fig = plt.gcf()
    plt.title(model + ' ' + polymer_name[0].upper() + polymer_name[1:] + ' ' + str_names.replace('_', ' '), fontsize=9)
    plt.scatter(x=results[:, 0], y=results[:, 1], cmap=plt.cm.get_cmap('nipy_spectral', len(polynucleotides[0])))
    plt.draw()
    plt.show()

    fig.savefig("data\\Cluster Analysis\\" + model + "_" + polymer_name + "_" + str_names + ".png", format='png')


def principal_component_analysis(df, genome_names, polymer_len):
    """Linear dimensionality reduction algorithm. Machine Learning algorithm for Cluster visualisation.
    Time in seconds, printed.

    :param DataFrame df: holds normalised frequency scores for fit & transformation
    :param list polynucleotides: combination of nucleic acids of n length each
    :param list genome_names: names extracted from included genomes['name']
    :param int polymer_name: len(polynucleotides[0]) used to get polymer name from POLYMER dict
    """
    time_start = time.time()
    pca = PCA(n_components=2)
    pca_results = pca.fit_transform(df)
    print("\nPCA complete. Time elapsed: {0:.2f}secs".format(time.time() - time_start))

    print("\n" + str(pca_results), type(pca_results))

    seaborn_scatterplot("PCA", pca_results, genome_names, polymer_len)


def tSNE(pca_results, genome_names, polymer_len):
    """t-distributed Stochastic Neighbour Embedding. Machine Learning algorithm for Cluster visualisation.
    Time in seconds, printed.

    :param ndarray pca_results: normalised frequency scores, reduced to 2 features, for fit & transformation
    :param list genome_names: names extracted from included genomes['name']
    :param int polymer_name: len(polynucleotides[0]) used to get polymer name from POLYMER dict
    """
    np.random.seed(42)  # Reproducability of results

    time_start = time.time()
    tsne = TSNE(n_components=2, perplexity=30, early_exaggeration=12)
    tsne_results = tsne.fit_transform(pca_results)
    print("\nt-SNE complete. Time elapsed: {0:.2f}secs".format(time.time() - time_start))

    print("\n" + str(tsne_results), type(tsne_results))

    seaborn_scatterplot("t-SNE", tsne_results, genome_names, polymer_len)
