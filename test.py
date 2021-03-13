import os
from os import listdir
from os.path import isfile, join
import random
import copy
import unittest
#from hypothesis import given, strategies as st
from threading import Thread
import Nucleotides
import Syntheses
import StandardFunctions as sf

BASES = ['A', 'C', 'G', 'T']
PROTEIN = ['I', 'M', 'T', 'N', 'K', 'S', 'R', 'L', 'P', 'H', 'Q', 'R', 'V', 'A', 'D', 'E', 'G', 'S', 'F', 'L', 'Y', 'C',
           'W', '*', '-']
DNA_SEQUENCE = ['A', 'G', 'T', 'C', 'C', 'A', 'G', 'T', 'G', 'T', 'A', 'A']
GENOME = {'name': "test", 'description': "test genome", 'sequence': DNA_SEQUENCE}
PREDICTED_OUTPUT = 'SPV*'


def file_exists(directory, filename):
    try:
        os.remove(directory + filename)
        flag = 1  # True - file did exist
    except FileNotFoundError:
        flag = 0  # False - file does not exist (Prediction)
    return flag


class TestNucleotides(unittest.TestCase):
    """ Nucleotides.py """

    def generate_sequence(bio_type):
        """Generates a pseudo-random DNA or Protein sequence for unbiased testing
        :param str bio_type: "dna" or "protein"

        :return: list test_sequence: pseudo-random sequence, based on average length for its 'bio_type'
        """
        print("1111111111111")
        if bio_type.lower() == "dna":
            items = BASES
            denominator = 2  # dimers
        elif bio_type.lower() == "protein":
            items = PROTEIN
            denominator = 1  # average protein sequence length
        else:
            print("Warning in test.py: biotype \"" + bio_type + "\" not recongnised. \nEnter \"DNA\" or \"Protein\"")
            return  # Exception Handling

        # Capture all filenames of 'bio_type'
        filenames = sf.capture_filenames(bio_type)
        directory = sf.directory(bio_type)

        lengths = []
        for f in filenames:
            with open(directory + f) as genome:
                content = sf.remove_firstline(directory + f)
                lengths.append(len(content.strip()))

        avg_len = sum(lengths) / len(lengths)  # test genome length
        num_polymers = int(avg_len / denominator)

        test_sequence = random.choices(items, k=num_polymers)
        print("Test Sequence:\n", test_sequence)

        return test_sequence

    def test_base_combinations_dimers(self):
        print("22222222222222222222")
        self.base_combinations = Nucleotides.base_combinations(2)  # dimers
        predicted_output = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG',
                            'TT']
        self.assertEqual(frozenset(self.base_combinations), frozenset(predicted_output),
                         'FAILED: Incorrect elemented generated for dimers list')  # frozenset() - alphabetical asc. order of lists

    def test_base_content(self, base_combinations):
        # Nucleotides.bases_content_plot("Test", )
        print("33333333333333333333333333")
        return

    #@given(st.floats(), st.floats())
    def test_normalised_frequencies(self, base_combinations, *genomes):
        """
        Genomes tend to have majority AT content over GC.

        """
        print("444444444444444444444444")
        # Nucleotides.composition_comparison(base_combinations, *genomes, test_genome)
        return

    def run_test_cases(self):
        #Thread(target=self.generate_sequence).start()
        Thread(target=self.test_base_combinations_dimers).start()
        #Thread(target=self.test_base_content).start()
        #Thread(target=self.test_normalised_frequencies).start()


class TestSyntheses(unittest.TestCase):
    """ Syntheses.py | protein() """

    def protein_sequence(self, name):
        file = open('data/Syntheses/' + 'protein_' + name, 'r')
        return file.read()

    def test_translation(self):
        # Conversion
        test_genome = copy.deepcopy(GENOME)
        Syntheses.protein(test_genome)
        ts = TestSyntheses()
        self.protein = ts.protein_sequence(test_genome['name'])
        self.assertEqual(self.protein, PREDICTED_OUTPUT, 'FAILED: Incorrect translation')

    def test_mononucleotide(self):
        # Non-trinculeotide composition (1 extra character)
        test_genome = copy.deepcopy(GENOME)
        test_genome['sequence'] += ['G']
        Syntheses.protein(test_genome)
        ts = TestSyntheses()
        self.protein = ts.protein_sequence(test_genome['name'])
        self.assertEqual(self.protein, PREDICTED_OUTPUT, 'FAILED: 1 additional DNA character not discarded')

    def test_dinucleotide(self):
        # Non-trinculeotide composition (2 extra characters)
        test_genome = copy.deepcopy(GENOME)
        test_genome['sequence'] += ['G', 'T']
        Syntheses.protein(test_genome)
        ts = TestSyntheses()
        self.protein = ts.protein_sequence(test_genome['name'])
        self.assertEqual(self.protein, PREDICTED_OUTPUT, 'FAILED: 2 additional DNA characters not discarded')

    def test_non_dna_alphabet(self):
        # Non-DNA alphabet character
        test_genome = copy.deepcopy(GENOME)
        test_genome['name'] = 'test_not_created_non_dna'
        test_genome['sequence'] += ['!']
        Syntheses.protein(test_genome)
        self.flag = file_exists('data/Syntheses/', 'protein_' + test_genome['name'])
        self.assertFalse(self.flag, 'FAILED: Protein Syntheses was not cancelled (non-DNA char)')

    def test_word(self):
        # Element longer than 1 character
        test_genome = copy.deepcopy(GENOME)
        test_genome['name'] = 'test_not_created_word'
        test_genome['sequence'] += ['WORD']
        Syntheses.protein(test_genome)
        self.flag = file_exists('data/Syntheses/', 'protein_' + test_genome['name'])
        self.assertFalse(self.flag, 'FAILED: Protein Syntheses was not cancelled (word)')

    def test_insufficient(self):
        test_genome = copy.deepcopy(GENOME)
        test_genome['name'] = 'test_not_created_insufficient'
        test_genome['sequence'] = test_genome['sequence'][:2]  # first 2 characters ['A', 'G']
        Syntheses.protein(test_genome)
        self.flag = file_exists('data/Syntheses/', 'protein_' + test_genome['name'])
        self.assertFalse(self.flag, 'FAILED: Protein Syntheses was not cancelled (insufficient)')

    def test_empty(self):
        test_genome = copy.deepcopy(GENOME)
        test_genome['name'] = 'test_not_created_empty'
        test_genome['sequence'] = []
        Syntheses.protein(test_genome)
        self.flag = file_exists('data/Syntheses/', 'protein_' + test_genome['name'])
        self.assertFalse(self.flag, 'FAILED: Protein Syntheses was not cancelled (empty)')

    def run_test_cases(self):
        Thread(target=self.test_translation).start()
        Thread(target=self.test_mononucleotide()).start()
        Thread(target=self.test_dinucleotide()).start()
        Thread(target=self.test_non_dna_alphabet).start()
        Thread(target=self.test_word).start()
        Thread(target=self.test_insufficient).start()
        Thread(target=self.test_empty).start()
        os.remove('data/Syntheses/protein_test')


class TestStandardFunctions(unittest.TestCase):
    """ StandardFunctions.py """

    def test_capture_filenames(self):
        predicted_dna_filenames = [d for d in listdir('src/') if isfile(join('src/', d))]
        predicted_protein_filenames = [p for p in listdir('data/Syntheses/') if isfile(join('data/Syntheses/', p))]
        self.dna_filenames = sf.capture_filenames("DNA")
        self.protein_filenames = sf.capture_filenames("Protein")
        self.assertFalse(sf.capture_filenames("foo"))  # no return
        self.assertNotEqual(self.dna_filenames, self.protein_filenames,
                            "FAILED: Same directory and files as each other")
        self.assertEquals(frozenset(predicted_dna_filenames), frozenset(self.dna_filenames),
                          "FAILED: wrong filenames (DNA)")  # frozenset() - alphabetical asc. order in dictionary
        self.assertEquals(frozenset(predicted_protein_filenames), frozenset(self.protein_filenames),
                          "FAILED: wrong filenames (Protein)")

    def test_output_name(self):
        test_genome_x = copy.deepcopy(GENOME)
        test_genome_x['name'] = 'foo'
        test_genome_y = copy.deepcopy(GENOME)
        test_genome_y['name'] = 'bar'
        test_genomes = (test_genome_x, test_genome_y)  # mimics tuples created by *args
        self.name = sf.output_name(test_genomes)
        self.assertEquals(self.name, "_foo_bar", "FAILED: incorrect output name returned")

    def test_output_filename(self):
        predicted_dna_filenames = [d for d in listdir('src/') if isfile(join('src/', d))]
        predicted_dna_filenames = [d.replace('.fasta', '') for d in predicted_dna_filenames]  # removes file extensions
        predicted_dna_filenames = [d.replace('.fna', '') for d in predicted_dna_filenames]
        predicted_protein_filenames = [p for p in listdir('data/Syntheses/') if isfile(join('data/Syntheses/', p))]
        predicted_protein_filenames = [d.replace('protein_', '') for d in predicted_protein_filenames]  # removes prefix
        file = open('data/Alignments/test.aln', 'w')
        # Capture all filenames
        bio_types = ["DNA", "Protein"]
        directories = ['src/', 'data/Syntheses/']
        output_filenames = []  # 2 lists - 'self.dna_filenames' & 'self.protein_filenames'
        for b, d, in zip(bio_types, directories):
            filenames = [f for f in listdir(d) if isfile(join(d, f))]
            output_filenames_bio_type = []
            for f in filenames:
                with open(d + f, 'r') as file:
                    name = sf.output_filename(f)  # Genome name
                output_filenames_bio_type.append(name)
            output_filenames.append(output_filenames_bio_type)  # append 'bio_type' list of 'filenames'
        file.close()
        os.remove('data/Alignments/test.aln')
        self.dna_filenames = output_filenames[0]
        self.protein_filenames = output_filenames[1]
        self.assertFalse(sf.capture_filenames("foo"))  # no return
        self.assertEquals(frozenset(self.dna_filenames), frozenset(predicted_dna_filenames),
                          "FAILED: wrong filenames (DNA)")  # frozenset() - alphabetical asc. order in dictionary
        self.assertEquals(self.protein_filenames, predicted_protein_filenames,
                          "FAILED: wrong filenames (Protein)")

    def test_ignore_firstline(self):
        test = open('test.fasta', 'w')
        test.write("first line" + "\n")
        test.write("second line" + "\n")
        test.close()
        self.content = sf.ignore_firstline('test.fasta')
        os.remove('test.fasta')
        self.assertNotEqual(self.content, "first linesecond line", "FAILED: first line not ignored")

    def test_directory(self):
        self.dir_dna = sf.directory("DNA")
        self.dir_protein = sf.directory("Protein")
        self.assertFalse(sf.directory("foo"))  # no return
        self.assertTrue(self.dir_dna, 'src/')
        self.assertTrue(self.dir_protein, 'data/Syntheses/')
        self.assertNotEqual(self.dir_dna, self.dir_protein, "FAILED: 'bio_type' directories shouldn't equal each other")

    def test_remove_prefix(self):
        self.name = sf.remove_prefix("foo_bar", "foo_")
        self.assertFalse(sf.remove_prefix("foo_", "foo_"))  # no return
        self.assertEqual(self.name, "bar", "FAILED: prefix not removed correctly")

    def run_test_cases(self):
        Thread(target=self.test_capture_filenames).start()
        Thread(target=self.test_output_name).start()
        Thread(target=self.test_output_filename).start()
        Thread(target=self.test_ignore_firstline).start()
        Thread(target=self.test_directory).start()
        Thread(target=self.test_remove_prefix).start()


def run():
    """ Instantiates Test Classes, run all methods """
    tn = TestNucleotides()
    tn.run_test_cases()

    # ts = TestSyntheses()
    # ts.run_test_cases()

    # tsf = TestStandardFunctions()
    # tsf.run_test_cases()
