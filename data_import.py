from os import listdir
from os.path import isfile, join


def read_genomes():
    """Imports Genome datasets from 'src/' folder.

    :return list: Stores complete genome datasets, each as dict
    """
    # Capture filenames from 'src/'
    filenames = []
    [filenames.append(g) for g in listdir('src/') if isfile(join('src/', g))]

    genomes = []  # 'name', 'description', 'sequence'

    for f in filenames:
        with open('src/' + f) as genome:
            description = genome.readline()  # '>'
            content = genome.read()

        name = f.rsplit('.', 1)[0]  # removes file extensions (eg. .fasta, .fna)
        name = name.split('.', 1)[0]  # removes gene version number of virus (eg. ".1")
        description = description.strip()  # removes newlines
        sequence = list(filter(lambda x: x != '\n', content))  # removes newlines

        g = {'name': name, 'description': description, 'sequence': sequence}
        genomes.append(g)

    print("Genomes: " + str(len(genomes)))

    return genomes