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

# Constants - appropriate referencing (print, charts, filenames)
polymer_type = {1: "monomers", 2: "dimers", 3: "trimers", 4: "tetramers", 5: "pentamers"}
path = 'data/Normalised Frequency/'
colormap = ["red", "green", "blue", "yellow", "purple"]
coronaviridae_names = []  # Appended to


def read_normalised_frequencies(polymers, name1, name2, *others):
    # Read in Normalised Frequency files (depending on: polymers & genomes of interest)
    # Minimum of 2 Coronaviridae names required - *others is optional
    coronaviridae_names.extend([name1, name2, *others])
    n = len(polymers)
    # "nf_" + polymer_type[n] + coronaviridae_names[name] + ".csv"

    # Capture filenames from data/Normalised Frequency
    # Append each file as a column
    nf_df = pd.DataFrame()  # DataFrame passed to PCA (dimentionality reduction algorithm)
    for file, cv_name in zip(os.listdir(path), coronaviridae_names):
        if os.path.isfile(os.path.join(path, file)):
            column = pd.read_csv(os.path.join(path, file))
            print(column)
            nf_df = pd.concat([nf_df, column], axis=1)

    print(nf_df)

    # Normalised Frequencies dataframe is passed through Dimentionality Reductionn algorithms for Cluster Analysis
    # Measuring time elapsed
    principal_component_analysis(nf_df)
    tSNE(nf_df)


def seaborn_scatterplot(model, results):
    # Machine Learning models output results
    labels = []
    for name in coronaviridae_names:
        labels.append(name * 100)

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
    fig.savefig("data\\Cluster Analysis\\" + model + "_Coronaviridae.png", format='png')


def principal_component_analysis(df):
    time_start = time.time()
    pca = PCA(n_components=2)
    pca_results = pca.fit_transform(df)
    print('\nPCA complete. Time elapsed: {} secs'.format(time.time() - time_start))

    print("\n" + str(pca_results), type(pca_results))

    seaborn_scatterplot("PCA", pca_results)


def tSNE(df):
    np.random.seed(42)  # Reproducability of results

    time_start = time.time()
    tsne = TSNE(n_components=2, perplexity=30, early_exaggeration=12)
    tsne_results = tsne.fit_transform(df)
    print('\nt-SNE complete. Time elapsed: {} secs'.format(time.time() - time_start))

    # print("\n" + tsne_results, type(tsne_results))

    seaborn_scatterplot("t-SNE", tsne_results)
