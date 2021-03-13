from Bio import Entrez
from urllib.request import HTTPError
import time
import re
from os import listdir
from os.path import isfile, join


def api(db, ids):
    """Invokes API, conditionally. Additionally, imports Genome datasets from 'src/' folder.

    :param str db: database name for type of data to extract
    :param list ids: list of str; IDs of datasets to download

    :return list genomes: Stores complete genome datasets, each as dict, from both sources
    """
    Entrez.email = "dbell30@qub.ac.uk"  # compulsory - identification for NCBI
    assert Entrez.email is not None

    genomes = []
    for id in ids:
        try:
            handle = Entrez.efetch(db=db, id=id, retmode='xml')
        except HTTPError:
            time.sleep(5)  # request too fast for server to respond
            handle = Entrez.efetch(db=db, id=id, retmode='xml')

        # Extraction
        response = Entrez.read(handle)  # 'Bio.Entrez.Parser.ListElement'
        name = response[0]['GBSeq_primary-accession']
        description = response[0]['GBSeq_definition']
        sequence = response[0]['GBSeq_sequence'].upper()
        g = {'name': name, 'description': description, 'sequence': sequence}
        genomes.append(g)

        # Create file
        if db == 'protein':
            file = open('data\\Syntheses\\' + 'protein_' + name, 'w')
            content = sequence
        else:  # e.g. DNA
            file = open('src\\' + name + '.fasta', 'w')
            content = '>' + name + ' ' + description + '\n'
            content += re.sub("(.{70})", "\\1\n", sequence, 0, re.DOTALL)  # new line every 70 chars

        file.write(content)
        file.close()

        time.sleep(1)

    return genomes


def read_genomes(db, ids):
    """Invokes API, conditionally. Additionally, imports Genome datasets from 'src/' folder.

    :param str db: database name for type of data to exract
    :param list ids: list of str; IDs of datasets to download

    :return list genomes: Stores complete genome datasets, each as dict, from both sources
    """
    genomes = []  # 'name', 'description', 'sequence'

    if (db != "") and (ids != ""): genomes = api(db, ids)

    # Capture filenames from 'src/'
    filenames = []
    [filenames.append(g) for g in listdir('src/') if isfile(join('src/', g))]

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
