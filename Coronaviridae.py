# Search Filters: United Kingdom, Human Host, Complete genome sequences
# https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=corona
from os import listdir
from os.path import isfile, join

# Capture filenames from src
filenames = []
test = [filenames.append(g) for g in listdir('src/') if isfile(join('src/', g))]


def readInGenomes(filenames):
    for f in filenames:
        with open('src/' + f) as genome:
            description = genome.readline() # '>'
            content = genome.read()

        print(description)
        sequence = list(filter(lambda a: a != '\n', content))  # removes newlines
        print(sequence)
        print('\n')

readInGenomes(filenames)

