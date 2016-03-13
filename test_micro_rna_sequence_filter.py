# -*- coding: utf-8 -*-


import unittest
import ttk
import Tkinter
from micro_rna_sequence_filter import SequenceFinder


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

        # self.program.entry = ["pozycja", "nukleotyd"]
        # self.program.positions = ["brak"] + range(1, 26)
        # self.program.nukleotydy = ["brak", "A", "G", "C", "T", "U"]

        self.program.open_input_file(self.program.input_file_a)

        self.program.position1 = ttk.Combobox()
        self.program.position2 = ttk.Combobox()
        self.program.position3 = ttk.Combobox()
        self.program.nukleotyd1 = ttk.Combobox()
        self.program.nukleotyd2 = ttk.Combobox()
        self.program.nukleotyd3 = ttk.Combobox()

    def test_start_SequenceFinder(self):
        """Testy parametrów uruchamiania programu"""

        self.assertEqual(self.program.entry, ["pozycja", "nukleotyd"])
        self.assertEqual(self.program.positions, ["brak"] + range(1, 26))
        self.assertEqual(self.program.nukleotydy, ["brak", "A", "G", "C", "T", "U"])

    def test_input_control(self):

        """ Testowanie danych input, parametrów przeszukiwania sekwencji."""
        self.program.position1.insert(0, self.program.entry[0])
        self.program.position2.insert(0, self.program.positions[0])
        self.program.position3.insert(0, 3)
        self.program.nukleotyd1.insert(0, self.program.entry[1])
        self.program.nukleotyd2.insert(0, self.program.nukleotydy[0])
        self.program.nukleotyd3.insert(0, "U")

        self.program.input_control()

        self.assertEquals(self.program.pos1, False,
                          msg="test_input_control - {}".format(self.program.pos1))
        self.assertEquals(self.program.nuc1, False,
                          msg="test_input_control - {}".format(self.program.nuc1))

        self.assertEquals(self.program.pos2, False,
                          msg="test_input_control - {}".format(self.program.pos2))
        self.assertEquals(self.program.nuc2, False,
                          msg="test_input_control - {}".format(self.program.nuc2))

        self.assertEquals(self.program.pos3, "3",
                          msg="test_input_control - {}".format(self.program.pos3))
        self.assertEquals(self.program.nuc3, "U",
                          msg="test_input_control - {}".format(self.program.nuc3))


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
        # sprawdzane są w poprzednich testach

        self.program.selekcja_sekwencji_1(0, "A")
        self.assertEqual(self.program.output_seq, expected_output_selekcja_1,
                         msg="Sprawdź --> test_selekcja_sekwencji_1")

        self.program.open_input_file(self.program.input_file_a)
        self.program.selekcja_sekwencji_2(0, "A", 1, "U")
        self.assertEqual(self.program.output_seq, expected_output_selekcja_2,
                         msg="Sprawdź --> test_selekcja_sekwencji_2")

        self.program.open_input_file(self.program.input_file_a)
        self.program.selekcja_sekwencji_3(0, "A", 1, "U", 2, "G")
        self.assertEqual(self.program.output_seq, expected_output_selekcja_3,
                         msg="Sprawdź --> test_selekcja_sekwencji_3")

    def test_go(self):
        """Test wyszukiwania sekwencji w przypadku fałaszywych danych wyszukiwania."""

        self.program.position1['values'] = "pozycja"
        self.program.position2['values'] = "pozycja"
        self.program.position3['values'] = "pozycja"
        self.program.nukleotyd1['values'] = "nukleotyd"
        self.program.nukleotyd2['values'] = "nukleotyd"
        self.program.nukleotyd3['values'] = "nukleotyd"
        self.program.text_field = Tkinter.Text()

        self.program.go()
        self.assertEqual(self.program.output_seq, "brak sekwencji spełniających zadane warunki")

if __name__ == '__main__':
    unittest.main()
