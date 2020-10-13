# Search Filters: United Kingdom, Human Host, Complete genome sequences
# https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=corona
from os import listdir
from os.path import isfile, join

# Capture filenames from src
filenames = []
test = [filenames.append(g) for g in listdir('src/') if isfile(join('src/', g))]


def readInGenomes(filenames):
    coronaviridae = [] # Virus name, description, complete genome sequence
    for f in filenames:
        with open('src/' + f) as genome:
            description = genome.readline() # '>'
            content = genome.read()

        coronavirus = f.rsplit('.', 1)[0]  # removes file extensions (eg. .fasta, .fna)
        description = description.strip() # removes newlines
        sequence = list(filter(lambda x: x != '\n', content))  # removes newlines

        coronavirius = { 'coronavirus': coronavirus, 'description': description, 'sequence': sequence }
        coronaviridae.append(coronavirius)

    print(coronaviridae)
    print("No. Coronaviridae Collected: " + str(len(coronaviridae)))
    return coronaviridae

readInGenomes(filenames)

