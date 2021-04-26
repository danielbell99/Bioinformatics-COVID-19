import numpy as np
import os
from os import listdir
from os.path import isfile, join
import random
import copy
import itertools
import time
from prettytable import PrettyTable
import unittest
from threading import Thread
import Nucleotides
import Syntheses
import PairwiseSequencing
from PairwiseSequenceAlignment import NeedlemanWunsch, SmithWaterman, tasks
import StandardFunctions as sf

BASES = ['A', 'C', 'G', 'T']
PROTEIN = ['I', 'M', 'T', 'N', 'K', 'S', 'R', 'L', 'P', 'H', 'Q', 'R', 'V', 'A', 'D', 'E', 'G', 'S', 'F', 'L', 'Y', 'C',
           'W', '*', '-']

DNA_SEQUENCE = ['A', 'G', 'T', 'C', 'C', 'A', 'G', 'T', 'G', 'T', 'A', 'A']
GENOME = {'name': "test", 'description': "test genome", 'sequence': DNA_SEQUENCE}
PREDICTED_OUTPUT = 'SPV*'

POINTS_SCHEME = {'match': 2.0, 'mismatch': -1.0, 'gap': -5.0}


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
        self.content = sf.polynucleotides(2)
        predicted = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
        self.assertEqual(frozenset(self.content), frozenset(predicted),
                         'FAILED: Incorrect elemented generated for dimers list')  # frozenset() - alphabetical asc. order of lists

    def test_normalised_frequencies_dimers(self):
        polynucleotides = sf.polynucleotides(2)  # dimer base combinations
        test_genome = copy.deepcopy(GENOME)
        test_sequence = ''.join(test_genome['sequence'])
        # Nomalised Frequencies - of their own dimer combination
        predicted_output = [0.17, 0.0, 0.33, 0.0, 0.17, 0.17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.17, 0.17, 0.17, 0.0]
        output = [round(nf, 2) for nf in Nucleotides.normalised_frequencies(polynucleotides, 'test', test_sequence)]
        os.remove('data/Normalised Frequency/dimers/nf_dimers_test.csv')
        self.assertEqual(output, predicted_output, 'FAILED: incorrect normalised frequencies')

    def run_test_cases(self):
        Thread(target=self.generate_sequence).start()
        Thread(target=self.test_base_combinations_dimers).start()
        Thread(target=self.test_normalised_frequencies_dimers).start()


class TestSyntheses(unittest.TestCase):
    """ Syntheses.py | protein() """

    def test_protein_sequence(self, name):
        file = open('data/Syntheses/' + 'protein_' + name, 'r')
        return file.read()

    def test_translation(self):
        test = open('src/' + 'test_translation' + '.fasta', 'w')
        test.write('description\n')
        test.write(''.join(DNA_SEQUENCE))
        test.close()
        # Conversion
        Syntheses.protein('test_translation')
        self.protein = sf.protein_sequence('test_translation')
        os.remove('src/' + 'test_translation.fasta')
        os.remove('data/Syntheses/protein_' + 'test_translation')
        self.assertEqual(self.protein, PREDICTED_OUTPUT, 'FAILED: Incorrect translation')

    def test_mononucleotide(self):
        test_genome = copy.deepcopy(GENOME)
        # File
        test = open('src/' + str(test_genome['name']) + '.fasta', 'w')
        test.write('decription\n')
        test.write(''.join(test_genome['sequence'] + ['G']))  # subject sequence
        test.close()
        # Non-trinculeotide composition (1 extra character)
        Syntheses.protein(str(test_genome['name']))
        self.protein = sf.protein_sequence(str(test_genome['name']))
        os.remove('src/' + str(test_genome['name'] + '.fasta'))
        os.remove('data/Syntheses/protein_' + str(test_genome['name']))
        self.assertEqual(self.protein, PREDICTED_OUTPUT, 'FAILED: 1 additional DNA character not discarded')

    def test_non_dna_alphabet(self):
        test_genome = copy.deepcopy(GENOME)
        # File
        test = open('src/' + str(test_genome['name']) + '.fasta', 'w')
        test.write('decription\n')
        test.write(''.join(test_genome['sequence'] + ['!']))  # subject sequence
        test.close()
        # Non-DNA alphabet character
        Syntheses.protein(str(test_genome['name']))
        self.protein = sf.protein_sequence(str(test_genome['name']))
        os.remove('src/' + str(test_genome['name'] + '.fasta'))
        os.remove('data/Syntheses/protein_' + str(test_genome['name']))
        self.flag = file_exists('data/Syntheses/', 'protein_' + test_genome['name'])
        self.assertFalse(self.flag, 'FAILED: Protein Syntheses was not cancelled (non-DNA char)')

    def test_word(self):
        test_genome = copy.deepcopy(GENOME)
        test_genome['name'] = "test_word"
        # File
        test = open('src/' + str(test_genome['name']) + '.fasta', 'w')
        test.write('decription\n')
        test.write(''.join(test_genome['sequence'] + ['WORD']))  # subject sequence
        test.close()
        # Element > 1 character
        Syntheses.protein(str(test_genome['name']))
        self.protein = sf.protein_sequence(str(test_genome['name']))
        os.remove('src/' + str(test_genome['name'] + '.fasta'))
        os.remove('data/Syntheses/protein_' + str(test_genome['name']))
        self.flag = file_exists('data/Syntheses/', 'protein_' + test_genome['name'])
        self.assertFalse(self.flag, 'FAILED: Protein Syntheses was not cancelled (word)')

    def test_insufficient(self):
        test_genome = copy.deepcopy(GENOME)
        # File
        test = open('src/' + str(test_genome['name']) + '.fasta', 'w')
        test.write('decription\n')
        test.write(''.join(test_genome['sequence'][:2]))  # subject sequence (first 2 characters ['A', 'G'])
        test.close()
        # Sequence less than 3 characters
        Syntheses.protein(str(test_genome['name']))
        self.protein = sf.protein_sequence(str(test_genome['name']))
        os.remove('data/Syntheses/protein_' + str(test_genome['name']))
        self.flag = file_exists('data/Syntheses/', 'protein_' + test_genome['name'])
        self.assertFalse(self.flag, 'FAILED: Protein Syntheses was not cancelled (insufficient)')

    def test_empty(self):
        test_genome = copy.deepcopy(GENOME)
        test_genome['name'] = 'test_not_created_empty'
        test_genome['sequence'] = []
        # File
        test = open('src/' + str(test_genome['name']) + '.fasta', 'w')
        test.write('decription\n')
        test.write(''.join(test_genome['sequence']))  # subject sequence
        test.close()
        # Empty
        Syntheses.protein(str(test_genome['name']))
        os.remove('src/' + str(test_genome['name'] + '.fasta'))
        self.flag = file_exists('data/Syntheses/', 'protein_' + test_genome['name'] + '.fasta')
        self.assertFalse(self.flag, 'FAILED: Protein Syntheses was not cancelled (empty)')

    def run_test_cases(self):
        Thread(target=self.test_translation).start()
        Thread(target=self.test_mononucleotide()).start()
        Thread(target=self.test_non_dna_alphabet).start()
        Thread(target=self.test_word).start()
        Thread(target=self.test_insufficient).start()
        Thread(target=self.test_empty).start()


class TestStandardFunctions(unittest.TestCase):
    """ StandardFunctions.py """

    def test_capture_filenames(self):
        predicted_dna_filenames = [d for d in listdir('src/') if isfile(join('src/', d)) and "test" not in d]
        predicted_protein_filenames = [p for p in listdir('data/Syntheses/') if isfile(join('data/Syntheses/', p)) and "test" not in p]
        self.dna_filenames = sf.capture_filenames("DNA")
        self.dna_filenames = [f for f in self.dna_filenames if "test" not in f]
        self.protein_filenames = sf.capture_filenames("Protein")
        self.protein_filenames = [f for f in self.protein_filenames if "test" not in f]
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
        test_genomes = (str(test_genome_x['name']), str(test_genome_y['name']))  # mimics tuples created by *args
        self.name = sf.output_name(test_genomes)
        self.assertEquals(self.name, "_foo_bar", "FAILED: incorrect output name returned")

    def test_output_filename(self):
        predicted_dna_filenames = [d for d in listdir('src/') if isfile(join('src/', d)) and "test" not in d]
        predicted_dna_filenames = [d.replace('.fasta', '') for d in predicted_dna_filenames]  # removes file extensions
        predicted_dna_filenames = [d.replace('.fna', '') for d in predicted_dna_filenames]
        predicted_protein_filenames = [p for p in listdir('data/Syntheses/') if isfile(join('data/Syntheses/', p)) and "test" not in p]
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
                if "test" not in name: output_filenames_bio_type.append(name)
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

    def test_dna_sequence(self):
        test = open('src/' + 'test_dna' + '.fasta', 'w')
        test.write("first line" + "\n")
        test.write("second line" + "\n")
        test.close()
        self.content = sf.dna_sequence('test_dna')
        os.remove('src/' + 'test_dna' + '.fasta')
        self.assertNotEqual(self.content, "first linesecond line", "FAILED: first line not ignored")

    def test_protein_sequence(self):
        self.directory_prefix = 'data/Syntheses/protein_'
        self.test_name_ps = 'test_protein_seq'
        test = open(self.directory_prefix + self.test_name_ps, 'w')
        test.write(PREDICTED_OUTPUT)
        test.close()
        self.content = sf.protein_sequence(self.test_name_ps)
        os.remove(self.directory_prefix + self.test_name_ps)
        self.assertEqual(str(self.content), str(PREDICTED_OUTPUT),
                         "FAILED: either incorrect protein sequence stored or collected")

    def test_directory(self):
        self.dir_dna = sf.directory("DNA")
        self.dir_protein = sf.directory("Protein")
        self.assertFalse(sf.directory("foo"))  # no return
        self.assertTrue(self.dir_dna, 'src/')
        self.assertTrue(self.dir_protein, 'data/Syntheses/')
        self.assertNotEqual(self.dir_dna, self.dir_protein, "FAILED: 'bio_type' directories shouldn't equal each other")

    def test_polynucleotides_dimers(self):
        self.dimers = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
        self.content = sf.polynucleotides(2)
        self.assertEqual(self.content, self.dimers, "FAILED: Polynucleotides not correctly generated")

    def test_remove_prefix(self):
        self.name = sf.remove_prefix("foo_bar", "foo_")
        self.assertFalse(sf.remove_prefix("foo_", "foo_"))  # no return
        self.assertEqual(self.name, "bar", "FAILED: prefix not removed correctly")

    def run_test_cases(self):
        Thread(target=self.test_capture_filenames).start()
        Thread(target=self.test_output_name).start()
        Thread(target=self.test_output_filename).start()
        Thread(target=self.test_dna_sequence).start()
        Thread(target=self.test_directory).start()
        Thread(target=self.test_polynucleotides_dimers).start()
        Thread(target=self.test_remove_prefix).start()
        Thread(target=self.test_protein_sequence).start()


class TestPairwiseSequenceAlignment(unittest.TestCase):
    """ NeedlemanWunsh.py, SmithWaterman.py, PairwiseSequencing.py, tasks.py """
    bio_types = ['DNA', 'Protein']
    genomes = ["SARS-CoV-JQ316196", "MERS-MT387202", "SARS-CoV-2-MT873892"]
    pairs = [list(i) for i in itertools.combinations(genomes, 2)]  # Generate unique pairs

    @staticmethod
    def time_algorithm(psa, bio_type, genomes):
        tic = time.process_time()
        psa(bio_type, genomes, **POINTS_SCHEME)  # PSA Algorithm invocation
        toc = time.process_time()
        return str('{:.4f}s'.format(toc - tic))

    def test_sequence_identity(self):
        seq_id, gapless_seq_id = 54.55, 66.67
        seq_identity = tasks.sequence_identity('ATATAAGGT-T', 'TTGT-AGATCT')
        self.assertEqual(round(float(seq_identity[0]), 2), seq_id, "FAILED: seq_id")
        self.assertEqual(round(float(seq_identity[1]), 2), gapless_seq_id, "FAILED: gapless_seq_id")
        self.assertNotEqual(seq_identity[0], seq_identity[1], "FAILED: possibly same variable passed back twice")

    def test_empy_matrix(self):
        m, n = 4, 4
        matrix = np.array(tasks.empty_matrix(m, n))  # 'shape'
        self.assertEqual(matrix.shape[0], m, "FAILED: wrong x-shape returned")
        self.assertEqual(matrix.shape[1], n, "FAILED: wrong y-shape returned")
        self.assertIsInstance(matrix, np.ndarray, "FAILED: empty matrix must be an ndarray")

    def test_match_points(self):
        points_scheme = list(POINTS_SCHEME.values())
        seqA, seqB = 'A', 'A'
        points = tasks.match(seqA, seqB, points_scheme)
        self.assertEqual(points, points_scheme[0], "FAILED: points incorrectly awarded from scheme")

    def test_mismatch_points(self):
        points_scheme = list(POINTS_SCHEME.values())
        seqA, seqB = 'A', 'B'
        points = tasks.match(seqA, seqB, points_scheme)
        self.assertEqual(points, points_scheme[1], "FAILED: points incorrectly awarded from scheme")

    def test_gap_points(self):
        points_scheme = list(POINTS_SCHEME.values())
        seqA, seqB = 'A', '-'
        points = tasks.match(seqA, seqB, points_scheme)
        self.assertEqual(points, points_scheme[1], "FAILED: points incorrectly awarded from scheme")

    def test_gap_function(self):
        penalty = tasks.gap_function(42, 0)
        self.assertEqual(penalty, 0, "FAILED: points incorrectly deducted")
        penalty = tasks.gap_function(42, 1)
        self.assertEqual(penalty, -1, "FAILED: points incorrectly deducted")
        penalty = tasks.gap_function(42, 2)
        self.assertEqual(penalty, -2.85, "FAILED: points incorrectly deducted")

    def test_score_alignment(self):
        align1, signs, align2, score = tasks.score_alignment('T-TGGAATATA', 'TCTAGA-TGTT', list(POINTS_SCHEME.values()))
        self.assertEqual(align1, 'ATATAAGGT-T', "FAILED: aligned sequence should be in reverse order")
        self.assertEqual(signs, '.|.| ||.| |', "FAILED: incorrect alignment signs")
        self.assertEqual(align2, 'TTGT-AGATCT', "FAILED: aligned should be in reverse order")
        self.assertEqual(score, '  Score=-1.0', "FAILED: wrong alignment score")

    def test_preparation(self):
        genomes = ["SARS-CoV-JQ316196", "SARS-CoV-2-MT873892"]
        seqA, seqB, match_points, mismatch_points, gap_points = tasks.preparation('DNA', genomes, POINTS_SCHEME)
        self.assertEqual(seqA, sf.dna_sequence(genomes[0])[:len(seqA)], "FAILED: DNA sequence")  # '[:len(seqA)]' - same length
        self.assertEqual(seqB, sf.dna_sequence(genomes[1])[:len(seqB)], "FAILED: DNA sequence")
        self.assertIsNotNone(match_points, "FAILED: no 'match_points' in 'points_scheme'")
        self.assertIsNotNone(mismatch_points, "FAILED: no 'mis_match' points in 'points_scheme'")
        self.assertIsNotNone(mismatch_points, "FAILED: no 'mis_match' points in 'points_scheme'")
        self.assertNotEqual(match_points, mismatch_points, "FAILED: 'match_points' == 'mismatch_points'")
        self.assertNotEqual(gap_points, mismatch_points, "FAILED: 'gap_points' == 'mismatch_points'")
        self.assertNotEqual(match_points, gap_points, "FAILED: 'match_points' == 'gap_points'")

    def time_all(self):
        print("-- Timing Pairwise Sequencing Algorithms --")
        for bio_type in self.bio_types:
            for genomes in self.pairs:  # list of lists
                nw_time = self.time_algorithm(NeedlemanWunsch.NeedlemanWunsch, bio_type, genomes)
                sw_time = self.time_algorithm(SmithWaterman.SmithWaterman, bio_type, genomes)
                ps_time = self.time_algorithm(PairwiseSequencing.run, bio_type, genomes)

            # Pretty Table
            table = PrettyTable(["Algorithm", "Time"])
            table.add_row(["Needleman-Wunsch", nw_time])
            table.add_row(["Smith-Waterman", sw_time])
            table.add_row(["Bio.pairwise2", ps_time])
            table.title = bio_type
            print(table, flush=True)

    def run_test_cases(self):
        Thread(target=self.test_sequence_identity).start()
        Thread(target=self.test_empy_matrix).start()
        Thread(target=self.test_match_points).start()
        Thread(target=self.test_mismatch_points).start()
        Thread(target=self.test_gap_points).start()
        Thread(target=self.test_gap_function).start()
        Thread(target=self.test_score_alignment).start()
        Thread(target=self.test_preparation).start()
        Thread(target=self.time_all).start()


def run():
    """ Instantiates Test Classes, run all methods """
    tn = TestNucleotides()
    tn.run_test_cases()

    ts = TestSyntheses()
    ts.run_test_cases()

    tsf = TestStandardFunctions()
    tsf.run_test_cases()

    tpsa = TestPairwiseSequenceAlignment()
    tpsa.run_test_cases()
