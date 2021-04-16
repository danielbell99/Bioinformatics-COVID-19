import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import time
import seaborn as sns
import matplotlib.pyplot as plt

# Constants - appropriate referencing (e.g. iteration, print, file name, visual)
POLYMER = {1: "monomers", 2: "dimers", 3: "trimers", 4: "tetramers", 5: "pentamers",
           6: "hexamers"}  # lower-case for file names
PATH = 'data/Cluster Analysis/reads/'
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
    for file, cv_name in zip(os.listdir('data/Normalised Frequency/' + polymer_name + '/'), genome_names):
        file_dir = 'data/Normalised Frequency/' + polymer_name + '/' + 'nf_' + polymer_name + '_' + cv_name + '.csv'
        if os.path.isfile(os.path.join(file_dir)):
            column = pd.read_csv(os.path.join(file_dir), header=None)
            nf_df = pd.concat([nf_df, column], axis=1)
    # nf_df = nf_df.T  # transpose - rows = viruses, cols = normalised freqs.

    #principal_component_analysis(nf_df, genome_names, polymer_len)
    #tSNE(nf_df, genome_names, polymer_len)


def seaborn_scatterplot(model, results, labels, genome_names):
    """Standardised method for using Seaborn's Scatterplot to visualise Unsupervised ML results.

    :param str model: String literal name of model for saving figure as filename
    :param ndarray results: holds 100 reads reduced to 2 columns, per genome
    :param list labels: a label is assigned to a datapoint, appears as a key in the legend
    :param list genome_names: genome names of interest
    """
    str_names = '_'.join(genome_names)

    plt.title(model + ' ' + str_names.replace('_', ' '), fontsize=8)
    sns_plot = sns.scatterplot(x=results[:, 0], y=results[:, 1], hue=labels)
    if model.lower() == 'pca':
        plt.xlabel(model + ' 1')
        plt.ylabel(model + ' 2')
    plt.draw()
    plt.show()

    fig = sns_plot.get_figure()
    fig.savefig('data/Cluster Analysis/' + model + '/' + model + '_' + str_names + '.png')


def principal_component_analysis(data, labels, genome_names):
    """Linear dimensionality reduction algorithm. Machine Learning algorithm for Cluster visualisation.
    Time in seconds, printed.

    :param ndarray data: holds 100 reads across 32 columns, per genome
    :param list labels: a label is assigned to an instance, appears as a key in the legend
    :param list genome_names: genome names of interest
    """
    time_start = time.time()
    pca = PCA(n_components=2)
    pca_results = pca.fit_transform(data)
    print("\nPCA complete. Time elapsed: {0:.2f}s".format(time.time() - time_start))

    print("\n" + str(pca_results), type(pca_results))

    seaborn_scatterplot("PCA", pca_results, labels, genome_names)


def tsne(data, labels, genome_names):
    """t-distributed Stochastic Neighbour Embedding. Machine Learning algorithm for Cluster visualisation.
    Time in seconds, printed.

    :param ndarray data: holds 'n_reads' across 32 columns, per genome
    :param list labels: a label is assigned to an instance, appears as a key in the legend
    :param list genome_names: genome names of interest
    """
    np.random.seed(42)  # Reproducibility of results

    time_start = time.time()
    tsne = TSNE(n_components=2, perplexity=30, early_exaggeration=12)
    tsne_results = tsne.fit_transform(data)
    print("\nt-SNE complete. Time elapsed: {0:.2f}s".format(time.time() - time_start))

    print("\n" + str(tsne_results), type(tsne_results))

    seaborn_scatterplot("t-SNE", tsne_results, labels, genome_names)


def kmeans(data, labels, genome_names):
    """k-Means Clustering - no. clusters defined by no. genomes.
    Elbow Method - Experiment w/ different 'n_clusters' and plot results.
    Machine Learning algorithm for Cluster visualisation.
    Time in seconds, printed.

    :param ndarray data: holds 'n_reads' across 32 columns, per genome
    :param list labels: a label is assigned to an instance, appears as a key in the legend
    :param list genome_names: genome names of interest
    """
    str_names = '_'.join(genome_names)

    # Elbow Method
    wcss = []  # sum of distance-squared of each data point, in a given cluster, w/ respect to their centroid
    for n_clusters in range(2, 11):  # 'n_clusters' of 2 to 10
        km = KMeans(n_clusters=n_clusters, init='k-means++', max_iter=300, n_init=10)  # init='k-means++' - prevents "Random Initialisation Trap"
        km.fit(data)
        wcss.append(km.inertia_)  # sum of distance-squared, per cluster

    plt.plot(range(2, 11), wcss)
    plt.title('Elbow Method ' + str_names.replace('_', ' '), size=7)
    plt.xlabel('Clusters')
    plt.ylabel('WCSS')
    plt.savefig('data/Cluster Analysis/k-Means/ElbowMethod_' + str_names + '.png')
    plt.show()

    # Cluster Analysis
    time_start = time.time()
    km = KMeans(n_clusters=len(genome_names), init='k-means++', max_iter=300, n_init=10, random_state=0)
    km_results = km.fit_transform(data)
    print("\nk-Means complete. Time elapsed: {0:.2f}s".format(time.time() - time_start))

    print("\n" + str(km_results), type(km_results))

    seaborn_scatterplot("k-Means", km_results, labels, genome_names)

    # Centroid Clustering
    plt.title('k-Means Centroids ' + str_names.replace('_', ' '), size=7)
    plt.scatter(data[:, 0], data[:, 1])
    plt.scatter(km.cluster_centers_[:, 0], km.cluster_centers_[:, 1], s=100, c='red')
    plt.savefig('data/Cluster Analysis/k-Means/Centroids_' + str_names + '.png')
    plt.show()


def run(n_reads, *names):
    """Imports reads from genomes of interest.
    Minimum of 2 Genome names required - '*others' is optional.

    :param int n_reads: no. reads/ instances per genome (typically 100)
    :param ndarray *names: array of str genome names ('*' >=1 genomes)
    """
    genome_names = [name for name in names]

    data, labels = np.zeros(shape=(n_reads, 32)), []  # data & labels passed to PCA & t-SNE models
    for read, name in zip(os.listdir(PATH), genome_names):
        file_dir = PATH + name + '.txt'
        if os.path.isfile(os.path.join(file_dir)):
            reads = np.loadtxt(file_dir)
            data = np.append(data, reads, axis=0)
            labels += [name for x in range(n_reads)]  # match 'n_reads' w/ n 'labels'
    data = data[n_reads:, :]  # first 'n_reads' rows are empty (so as to establish a shape for appending)
    labels = np.array(labels)

    principal_component_analysis(data, labels, genome_names)
    tsne(data, labels, genome_names)
    kmeans(data, labels, genome_names)
