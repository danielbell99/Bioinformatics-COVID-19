import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import seaborn as sns
import matplotlib.pyplot as plt
import time

# Constants - appropriate referencing (e.g. iteration, print, file name, visual)
POLYMER = {1: "monomers", 2: "dimers", 3: "trimers", 4: "tetramers", 5: "pentamers"}  # lower-case for file names
PATH = 'data/Normalised Frequency/'
COLOUR_MAP = ["red", "green", "blue", "yellow", "purple"]


coronaviridae_names = []  # Appended to


def read_normalised_frequencies(polymers, name1, name2, *others):
    """Imports produced Normalised Frequency files (depending on: polymers & genomes of interest).
    Normalised Frequencies are passed through Dimentionality Reduction algorithms for Cluster Analysis.
    Measuring time elapsed.
    Minimum of 2 Coronaviridae names required - '*others' is optional.

    :param list polymers: combination of nucleic acids of n length each
    :param dict (name1, name2, *others): contains entire genome data; 'name' value used (*others is optional)
    """
    coronaviridae_names.extend([name1, name2, *others])
    n = len(polymers)
    # "nf_" + POLYMER[n] + coronaviridae_names[name] + ".csv"

    # Capture filenames from data/Normalised Frequency
    # Append each file as a column
    nf_df = pd.DataFrame()  # DataFrame passed to PCA (dimentionality reduction algorithm)
    for file, cv_name in zip(os.listdir(PATH), coronaviridae_names):
        if os.path.isfile(os.path.join(PATH, file)):
            column = pd.read_csv(os.path.join(PATH, file))
            print(column)
            nf_df = pd.concat([nf_df, column], axis=1)

    print(nf_df)

    principal_component_analysis(nf_df)
    tSNE(nf_df)


def seaborn_scatterplot(model, results):
    """Standardised method for using Seaborn's Scatterplot to visualise Machine Learning results.

    :param str model: String literal name of model for saving figure as filename
    :param ndarray results: contains entire genome data; 'name' value used (*others is optional)
    """
    # Machine Learning models output results
    labels = []
    [labels.append(name * 100) for name in coronaviridae_names]

    label_vec = ["SARS-CoV"] * 100
    label_vec = label_vec + ["bat-SL-CoV"] * 100
    label_vec = label_vec + ["COVID-19"] * 100

    fig = plt.figure(figsize=(7, 7))
    sns.scatterplot(x=results[:, 0],
                    y=results[:, 1],
                    # hue=label_vec,
                    alpha=1,
                    s=100,
                    palette=['red', 'green', 'blue']).plot()
    # fig = plt.gcf()
    # plt.scatter(x=results[:, 0], y=results[:, 1])
    # plt.show()
    # plt.draw()
    fig.savefig('data\\Cluster Analysis\\' + model + '_Coronaviridae.png', format='png')


def principal_component_analysis(df):
    """Linear dimensionality reduction algorithm. Machine Learning algorithm for Cluster visualisation.
    Time in seconds, printed.

    :param DataFrame df: holds normalised frequency scores for fit & transformation
    """
    time_start = time.time()
    pca = PCA(n_components=2)
    pca_results = pca.fit_transform(df)
    print("\nPCA complete. Time elapsed: {} secs".format(time.time() - time_start))

    print("\n" + str(pca_results), type(pca_results))

    seaborn_scatterplot("PCA", pca_results)


def tSNE(pca_results):
    """t-distributed Stochastic Neighbour Embedding. Machine Learning algorithm for Cluster visualisation.
    Time in seconds, printed.

    :param ndarray pca_results: normalised frequency scores, reduced to 2 features, for fit & transformation
    """
    np.random.seed(42)  # Reproducability of results

    time_start = time.time()
    tsne = TSNE(n_components=2, perplexity=30, early_exaggeration=12)
    tsne_results = tsne.fit_transform(pca_results)
    print("\nt-SNE complete. Time elapsed: {} secs".format(time.time() - time_start))

    # print("\n" + tsne_results, type(tsne_results))

    seaborn_scatterplot("t-SNE", tsne_results)
