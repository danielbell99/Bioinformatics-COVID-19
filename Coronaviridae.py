# Search Filters: United Kingdom, Human Host, Complete genome sequences
# https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=corona
from os import listdir
from os.path import isfile, join

def readInGenomes():
    # Capture filenames from src
    filenames = []
    [filenames.append(g) for g in listdir('src/') if isfile(join('src/', g))]

    coronaviridae = [] # Virus name, description, complete genome sequence

    for f in filenames:
        with open('src/' + f) as genome:
            description = genome.readline() # '>'
            content = genome.read()

        name = f.rsplit('.', 1)[0]  # removes file extensions (eg. .fasta, .fna)
        name = name.split('.', 1)[0]  # removes gene version number of virus eg. ".1"
        description = description.strip() # removes newlines
        sequence = list(filter(lambda x: x != '\n', content))  # removes newlines

        coronavirius = { 'name': name, 'description': description, 'sequence': sequence }
        coronaviridae.append(coronavirius)

    print("Coronaviridae Genomes: " + str(len(coronaviridae)))
    return coronaviridae
