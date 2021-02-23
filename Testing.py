import os
import random
import copy
import unittest
# from hypothesis import given, strategies as st
import threading
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
    # test_genome['name']
    try:
        # self.protein = ts.protein_sequence(test_genome['name'])
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
            print("Warning in Testing.py: biotype \"" + bio_type + "\" not recongnised. \nEnter \"DNA\" or \"Protein\"")
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
        self.base_combinations = Nucleotides.bases_combinations(2)  # dimers
        predicted_output = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG',
                            'TT']
        self.assertEqual(self.base_combinations, predicted_output,
                         'FAILED: Wrong Dimers generated, either ordering or different no. elements')

    def test_base_content(self, base_combinations):
        # Nucleotides.bases_content_plot("Test", )
        print("33333333333333333333333333")
        return

    # @given(st.floats(), st.floats())
    def test_normalised_frequencies(self, base_combinations, *genomes):
        """
        Coronaviridae tend to have majority AT content over GC.

        """
        print("444444444444444444444444")
        # Nucleotides.composition_comparison(base_combinations, *genomes, test_genome)
        return

    def runall(self):
        Thread(target=self.generate_sequence).start()
        Thread(target=self.test_base_combinations_dimers).start()
        Thread(target=self.test_base_content).start()
        Thread(target=self.test_normalised_frequencies).start()


class TestSyntheses(unittest.TestCase):
    """ Syntheses.py """

    def protein_sequence(self, name):
        file = open('data/Syntheses/' + 'protein_' + name, 'r')
        return file.read()

    def test_translation(self):
        print("5555555555555555555")
        # Conversion
        test_genome = copy.deepcopy(GENOME)
        Syntheses.protein(test_genome)
        ts = TestSyntheses()
        self.protein = ts.protein_sequence(test_genome['name'])
        self.assertEqual(self.protein, PREDICTED_OUTPUT, 'FAILED: Incorrect translation')

    def test_mononucleotide(self):
        print("6666666666666666")
        # Non-trinculeotide composition (1 extra character)
        test_genome = copy.deepcopy(GENOME)
        test_genome['sequence'] += ['G']
        Syntheses.protein(test_genome)
        ts = TestSyntheses()
        self.protein = ts.protein_sequence(test_genome['name'])
        self.assertEqual(self.protein, PREDICTED_OUTPUT, 'FAILED: 1 additional DNA character not discarded')

    def test_dinucleotide(self):
        print("777777777777777777")
        # Non-trinculeotide composition (2 extra characters or 1 less)
        test_genome = copy.deepcopy(GENOME)
        test_genome['sequence'] += ['G', 'T']
        print("#2", test_genome)
        Syntheses.protein(test_genome)
        ts = TestSyntheses()
        self.protein = ts.protein_sequence(test_genome['name'])
        self.assertEqual(self.protein, PREDICTED_OUTPUT, 'FAILED: 2 additional DNA characters not discarded')

    def test_non_dna_alphabet(self):
        print("888888888888888888")
        # Non-DNA alphabet character
        test_genome = copy.deepcopy(GENOME)
        test_genome['name'] = 'test_not_created_non_dna'
        test_genome['sequence'] += ['!']
        Syntheses.protein(test_genome)
        self.flag = file_exists('data/Syntheses/', 'protein_' + test_genome['name'])
        self.assertEqual(self.flag, 0, 'FAILED: Protein Syntheses was not cancelled (non-DNA char)')

    def test_word(self):
        print("9999999999999999999999")
        # Element longer than 1 character
        test_genome = copy.deepcopy(GENOME)
        test_genome['name'] = 'test_not_created_word'
        test_genome['sequence'] += ['WORD']
        Syntheses.protein(test_genome)
        self.flag = file_exists('data/Syntheses/', 'protein_' + test_genome['name'])
        self.assertEqual(self.flag, 0, 'FAILED: Protein Syntheses was not cancelled (word)')

    def test_insufficient(self):
        print("AAAAAAAAAAAAAAAAAAAA")
        test_genome = copy.deepcopy(GENOME)
        test_genome['name'] = 'test_not_created_insufficient'
        test_genome['sequence'] = test_genome['sequence'][:2]  # first 2 characters ['A', 'G']
        Syntheses.protein(test_genome)
        self.flag = file_exists('data/Syntheses/', 'protein_' + test_genome['name'])
        self.assertEqual(self.flag, 0, 'FAILED: Protein Syntheses was not cancelled (insufficient)')

    def test_empty(self):
        print("BBBBBBBBBBBBBBBBBBBBBBBB")
        test_genome = copy.deepcopy(GENOME)
        test_genome['name'] = 'test_not_created_empty'
        test_genome['sequence'] = []
        Syntheses.protein(test_genome)
        self.flag = file_exists('data/Syntheses/', 'protein_' + test_genome['name'])
        self.assertEqual(self.flag, 0, 'FAILED: Protein Syntheses was not cancelled (empty)')

    def runall(self):
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

    def runall(self):
        Thread(target=self.method).start()
        Thread(target=self.method).start()
        Thread(target=self.method).start()
        Thread(target=self.method).start()
        Thread(target=self.method).start()
        Thread(target=self.method).start()
        Thread(target=self.method).start()
        Thread(target=self.method).start()


def run():
    """ Instantiates Test Classes, run all methods """
    # tn = TestNucleotides()
    # tn.runall()

    ts = TestSyntheses()
    ts.runall()

    # tsf = TestStandardFunctions()
    # tsf.runall()
