# -*- coding: utf-8 -*-

import unittest
from MicroRNASequenceFilter import SequenceFinder

__author__ = 'Marcin Pieczyński'


###########################################
# Do not change text below. These strings were prepared only for testing !!!
###########################################
expected_false_output = """brak sekwencji spełniających zadane warunki"""

expected_output_selekcja_1 = """Dane z pliku test_file.txt
Work finished -> 8 sequences were found with:
A in position 1

stu-miR8008a
AUUUCCAGAAAAGCGACGGACAGU
stu-miR8005b-5p
ACUCUAAAUUUUAAAUUCUAAAUC
stu-miR1919-3p
ACGAGAGUCAUCUGUGACAGG
stu-miR8032a-3p
AGUGUGAGUCGGUGUGAUUAGG
stu-miR166d-5p
AGAAUGUCGUCUGGUUCGAGA
stu-miR1886b
AUGGUAUCGUGAGAUGAAAUCAGC
stu-miR7999-3p
ACGACCCGUAGAACUGCCCACGAC
stu-miR1886i-5p
AUGAGAUGAAAUUAGCGUUUGGAU
"""

expected_output_selekcja_2 = """Dane z pliku test_file.txt
Work finished -> 3 sequences were found with:
A in position 1
U in position 2

stu-miR8008a
AUUUCCAGAAAAGCGACGGACAGU
stu-miR1886b
AUGGUAUCGUGAGAUGAAAUCAGC
stu-miR1886i-5p
AUGAGAUGAAAUUAGCGUUUGGAU
"""


expected_output_selekcja_3 = """Dane z pliku test_file.txt
Work finished -> 2 sequences were found with:
A in position 1
U in position 2
G in position 3

stu-miR1886b
AUGGUAUCGUGAGAUGAAAUCAGC
stu-miR1886i-5p
AUGAGAUGAAAUUAGCGUUUGGAU
"""
############################################################


class MicroRNA_Sequence_Filter_Test(unittest.TestCase):
    """Tests for MicroRNA_Sequence_Filter2.py"""

    def setUp(self):
        """Tworzenie instancji dla klascy SequenceFinder"""
        self.program = SequenceFinder()
        self.program.input_file_a = "test_file.txt"
        self.program.open_input_file(self.program.input_file_a)

    def test_start_SequenceFinder(self):
        """Testy parametrów uruchamiania programu"""

        self.assertEqual(self.program.entry, ["pozycja", "nukleotyd"])
        self.assertEqual(self.program.positions, ["brak"] + range(1, 26))
        self.assertEqual(self.program.nukleotydy, ["brak", "A", "G", "C", "T", "U"])

    def test_input_control_1a(self):
        """Test 1a - testowanie parametrów przeszukiwania sekwencji - pierwszy nukleotyd"""

        self.program.input_control1("pozycja", "nukleotyd")
        self.assertEquals(self.program.pos1, False,
                          msg="test_input_control-1a - {}".format(self.program.pos1))
        self.assertEquals(self.program.nuc1, False,
                          msg="test_input_control-1a - {}".format(self.program.nuc1))

        self.program.input_control1("pozycja", "brak")
        self.assertEquals(self.program.pos1, False,
                          msg="test_input_control-1b - {}".format(self.program.pos1))
        self.assertEquals(self.program.nuc1, False,
                          msg="test_input_control-1b - {}".format(self.program.nuc1))

        self.program.input_control1("brak", "nukleotyd")
        self.assertEquals(self.program.pos1, False,
                          msg="test_input_control-1c - {}".format(self.program.pos1))
        self.assertEquals(self.program.nuc1, False,
                          msg="test_input_control-1c - {}".format(self.program.nuc1))

        self.program.input_control1("brak", "brak")
        self.assertEquals(self.program.pos1, False,
                          msg="test_input_control-1d - {}".format(self.program.pos1))
        self.assertEquals(self.program.nuc1, False,
                          msg="test_input_control-1d - {}".format(self.program.nuc1))

    def test_input_control_1b(self):
        """Test 1b - testowanie parametrów przeszukiwania sekwencji - pierwszy nukleotyd"""

        for x in self.program.positions[1:]:
            for y in self.program.nukleotydy[1:]:
                self.program.input_control1(x, y)
                self.assertEquals(self.program.pos1, x,
                                  msg="test_input_contro2 - {}-->{}".format(self.program.pos1, x))
                self.assertEquals(self.program.nuc1, y,
                                  msg="test_input_contro2 - {}-->{}".format(self.program.nuc1, y))

    def test_input_control_2a(self):
        """Test 2a - testowanie parametrów przeszukiwania sekwencji - drugi nukleotyd"""

        self.program.input_control2("pozycja", "nukleotyd")
        self.assertEquals(self.program.pos2, False,
                          msg="test_input_control-1a - {}".format(self.program.pos2))
        self.assertEquals(self.program.nuc2, False,
                          msg="test_input_control-1a - {}".format(self.program.nuc2))

        self.program.input_control1("pozycja", "brak")
        self.assertEquals(self.program.pos2, False,
                          msg="test_input_control-1b - {}".format(self.program.pos2))
        self.assertEquals(self.program.nuc2, False,
                          msg="test_input_control-1b - {}".format(self.program.nuc2))

        self.program.input_control1("brak", "nukleotyd")
        self.assertEquals(self.program.pos2, False,
                          msg="test_input_control-1c - {}".format(self.program.pos2))
        self.assertEquals(self.program.nuc2, False,
                          msg="test_input_control-1c - {}".format(self.program.nuc2))

        self.program.input_control1("brak", "brak")
        self.assertEquals(self.program.pos2, False,
                          msg="test_input_control-1d - {}".format(self.program.pos2))
        self.assertEquals(self.program.nuc2, False,
                          msg="test_input_control-1d - {}".format(self.program.nuc2))

    def test_input_control_2b(self):
        """Test 2b - testowanie parametrów przeszukiwania sekwencji - drugi nukleotyd"""

        for x in self.program.positions[1:]:
            for y in self.program.nukleotydy[1:]:
                self.program.input_control2(x, y)
                self.assertEquals(self.program.pos2, x,
                                  msg="test_input_contro2 - {}-->{}".format(self.program.pos2, x))
                self.assertEquals(self.program.nuc2, y,
                                  msg="test_input_contro2 - {}-->{}".format(self.program.nuc2, y))

    def test_input_control_3a(self):
        """Test 3a - testowanie parametrów przeszukiwania sekwencji - trzeci nukleotyd"""

        self.program.input_control3("pozycja", "nukleotyd")
        self.assertEquals(self.program.pos3, False,
                          msg="test_input_control-1a - {}".format(self.program.pos3))
        self.assertEquals(self.program.nuc3, False,
                          msg="test_input_control-1a - {}".format(self.program.nuc3))

        self.program.input_control1("pozycja", "brak")
        self.assertEquals(self.program.pos3, False,
                          msg="test_input_control-1b - {}".format(self.program.pos3))
        self.assertEquals(self.program.nuc3, False,
                          msg="test_input_control-1b - {}".format(self.program.nuc3))

        self.program.input_control1("brak", "nukleotyd")
        self.assertEquals(self.program.pos3, False,
                          msg="test_input_control-1c - {}".format(self.program.pos3))
        self.assertEquals(self.program.nuc3, False,
                          msg="test_input_control-1c - {}".format(self.program.nuc3))

        self.program.input_control1("brak", "brak")
        self.assertEquals(self.program.pos3, False,
                          msg="test_input_control-1d - {}".format(self.program.pos3))
        self.assertEquals(self.program.nuc3, False,
                          msg="test_input_control-1d - {}".format(self.program.nuc3))

    def test_input_control_3b(self):
        """Test 3b - testowanie parametrów przeszukiwania sekwencji - trzeci nukleotyd"""

        for x in self.program.positions[1:]:
            for y in self.program.nukleotydy[1:]:
                self.program.input_control3(x, y)
                self.assertEquals(self.program.pos3, x,
                                  msg="test_input_contro2 - {}-->{}".format(self.program.pos3, x))
                self.assertEquals(self.program.nuc3, y,
                                  msg="test_input_contro2 - {}-->{}".format(self.program.nuc3, y))

    def test_get_pos_from_str(self):
        """Test get_pos_from_str - testuje czy pozycja nukleotydu jest int i ma włąściwą wartość"""
        for x in self.program.positions[1:]:
            for y in self.program.positions[1:]:
                for z in self.program.positions[1:]:
                    self.program.pos1, self.program.pos2, self.program.pos3 = x, y, z
                    self.program.get_pos_from_str()
                    self.assertEquals(self.program.pos1, int(x)-1,
                                      msg="test_get_pos_from_str - {}-->{}".format(self.program.pos1, x-1))
                    self.assertEquals(self.program.pos2, int(y)-1,
                                      msg="test_get_pos_from_str - {}-->{}".format(self.program.pos2, y-1))
                    self.assertEquals(self.program.pos3, int(z)-1,
                                      msg="test_get_pos_from_str - {}-->{}".format(self.program.pos3, z-1))

    def test_check_brak_seq(self):
        """ Test sprawdzający brak wyników wyszukiwania"""
        self.program.output_seq = ""
        self.program.check_brak_seq()
        self.assertEqual(self.program.output_seq, "\nBrak sekwencji spełniających zadane warunki",
                         msg="Sprawdź --> test_check_brak_seq")

    def test_selekcja_sekwencji(self):
        """Test wyszukiwania sekwencji wedłóg zadanych parametrów dla i nukleotydu"""

        # Testowanie tylko prawidłowych danych dla self.program.pos i self.program.nuc
        # Przypadki z danymi nieprawidłowymi ("pozycja", "nukleotyd", "brak")
        # sprawdzane są pośrednio w poprzednich testach i
        # nie powinny dość do etapy wykonywania metody selekcja_sekwencji_1  / _2 / lub _3

        self.program.selekcja_sekwencji_1(0, "A")
        self.assertEqual(self.program.output_seq, expected_output_selekcja_1,)
                         # msg="Sprawdź --> test_selekcja_sekwencji_1")

        self.program.open_input_file(self.program.input_file_a)
        self.program.selekcja_sekwencji_2(0, "A", 1, "U")
        self.assertEqual(self.program.output_seq, expected_output_selekcja_2,
                         msg="Sprawdź --> test_selekcja_sekwencji_2")

        self.program.open_input_file(self.program.input_file_a)
        self.program.selekcja_sekwencji_3(0, "A", 1, "U", 2, "G")
        self.assertEqual(self.program.output_seq, expected_output_selekcja_3,
                         msg="Sprawdź --> test_selekcja_sekwencji_3")

if __name__ == '__main__':
    unittest.main()
