# -*- coding: utf-8 -*-

import Tkinter
import tkMessageBox
import tkFileDialog
import os
import ttk
from Bio import SeqIO

__author__ = 'Marcin Pieczyński'


"""
Program do przeszukiwania sekwencji fasta wedłóg zadanych parametrów.
Umożliwia selekcje sekwencji na podstawie wybranych nukleotydów w wybranych pozycjach.
"""

text_output1 = """Dane z pliku {}
Work finished -> {} sequences were found with:
{} in position {}\n\n"""

text_output2 = """Dane z pliku {}
Work finished -> {} sequences were found with:
{} in position {}
{} in position {}\n\n"""

text_output3 = """Dane z pliku {}
Work finished -> {} sequences were found with:
{} in position {}
{} in position {}
{} in position {}\n\n"""


class SequenceFinder(object):
    """
    Klasa zawierająca cały program, GUI oraz wszystkie metody.
    """
    def __init__(self):
        self.top = Tkinter.Tk()
        self.top.wm_title("MicroRNASequenceFilter")
        self.intro = "\nAplikacja do filtrowania sekwencji\n"
        self.entry = ["pozycja", "nukleotyd"]
        self.positions = ["brak"] + range(1, 26)
        self.nukleotydy = ["brak", "A", "G", "C", "T", "U"]
        self.filepath_input = Tkinter.StringVar()

        self.position1 = ttk.Combobox(self.top)
        self.position2 = ttk.Combobox(self.top)
        self.position3 = ttk.Combobox(self.top)
        self.nukleotyd1 = ttk.Combobox(self.top)
        self.nukleotyd2 = ttk.Combobox(self.top)
        self.nukleotyd3 = ttk.Combobox(self.top)

        self.output_seq = ""
        self.prefix = ""
        self.n = 0
        self.pos1 = None
        self.pos2 = None
        self.pos3 = None
        self.nuc1 = None
        self.nuc2 = None
        self.nuc3 = None
        self.input_file_a = None
        self.input_seq = None
        self.output_file_a = None
        self.output_file = None

    def open_input_file(self, in_file):
        """
        Funkcja otwierająca zadany plik z danymi wejściowymi- in_file.
        """
        try:
            self.input_seq = open(in_file, "rU")
        except IOError:
            return self.error_inputfile()

    def input_file(self):
        """
        Okno wprowadzające plik z danymi wejściowymi
        """
        try:
            f = tkFileDialog.askopenfilename(parent=self.top, initialdir="/home/",
                                             title="Wybór pliku fasta",
                                             filetypes=[("Text file", ".txt"), ("Fasta file", ".fasta")])
            self.filepath_input.set(os.path.realpath(f))
            self.input_file_a = os.path.realpath(f)
        except ValueError:
            return self.error_inputfilepath()

    def save(self):
        """
        Okno służące do zapisywania danych wyjściowych do pliku.
        """
        try:
            f = tkFileDialog.asksaveasfilename(parent=self.top, initialdir="/home/marcin/pulpit/Py/",
                                               title="Wybór pliku fasta", filetypes=[("Text file", ".txt")])
            self.output_file_a = os.path.realpath(f) + ".txt"
            self.output_file = open(self.output_file_a, "w")
            self.output_file.write(str(self.output_seq))
            self.output_file.close()
            self.work_finished()
        except ValueError:
            return self.error_outputfile()

    @staticmethod
    def error_inputfile():
        """
        Okno z informacją o błędzie podczas wczytywania pliku.
        """
        tkMessageBox.showerror("Error",
                               " Hej, jest problem z otwarciem pliku co go to chcesz załadować. Zrób z tym coś.")

    @staticmethod
    def error_inputfilepath():
        """
        Okno z informacją o nierawidłowej ścieżce do pliku.
        """
        tkMessageBox.showerror("Error",
                               " Hej, jest problem z ścieżką do pliku co go to chcesz załadować. Zrób z tym coś.")

    @staticmethod
    def error_outputfile():
        """
        Okno z informacją o błędzie podczas zapisywania pliku.
        """
        tkMessageBox.showerror("Error",
                               " Hej, jest problem z tym plikiem co go to chcesz zapisać. Zrób z tym coś.")

    @staticmethod
    def work_finished():
        """
        Okno z informacją o zapisaniu pliku.
        """
        tkMessageBox.showinfo("Wszystko OK", "Hej, wszystko OK, plik został zapisany we wskazanym miejscu :-)")

    def input_control(self):
        """
        Metoda sprawdzająca dane input dla pierwszego nukleotydu.
        """
        self.pos1, self.nuc1 = self.position1.get(), self.nukleotyd1.get()
        self.pos2, self.nuc2 = self.position2.get(), self.nukleotyd2.get()
        self.pos3, self.nuc3 = self.position3.get(), self.nukleotyd3.get()

        if self.pos1 == self.entry[0] or self.pos1 == self.positions[0] or \
           self.nuc1 == self.entry[1] or self.nuc1 == self.nukleotydy[0]:
            self.pos1, self.nuc1 = False, False

        if self.pos2 == self.entry[0] or self.pos2 == self.positions[0] or \
           self.nuc2 == self.entry[1] or self.nuc2 == self.nukleotydy[0]:
            self.pos2, self.nuc2 = False, False

        if self.pos3 == self.entry[0] or self.pos3 == self.positions[0] or \
           self.nuc3 == self.entry[1] or self.nuc3 == self.nukleotydy[0]:
            self.pos3, self.nuc3 = False, False

    def get_pos_from_str(self):
        """
        Metoda zamieniająca pozycje nukleotydów string --> int.
        """
        if self.pos1:
            self.pos1 = int(self.pos1)-1
        if self.pos2:
            self.pos2 = int(self.pos2)-1
        if self.pos3:
            self.pos3 = int(self.pos3)-1

    def check_brak_seq(self):
        """
        Metoda zwraca informacje o braku znalezionych sekwencji spełniających zadana warunki.
        """
        if not self.output_seq:
            self.output_seq = "\nBrak sekwencji spełniających zadane warunki"

    def start_output(self):
        """
        Metoda rozpoczynająca przeszukiwanie sekwencji, przypisywanie zmiennych na start.
        """
        self.output_seq, self.prefix, self.n = "", "", 0

    def selekcja(foo):
        """
        Metoda dokonująca właściwego filtrowania danych.
        """
        def wraper(self, *args):
            self.start_output()
            foo(self, *args)
            self.check_brak_seq()
            self.output_seq = self.prefix + self.output_seq
        return wraper

    @selekcja
    def selekcja_sekwencji_1(self, pozycja1, nukleotyd1):
        """
        Funkcja filtrująca sekwencje wedłóg pozycji jednego nukleotydu.
        """
        for record in SeqIO.parse(self.input_seq, "fasta"):
            if record.seq[pozycja1].upper() == nukleotyd1:
                self.output_seq += "{}\n{}\n".format(record.id, record.seq)
                self.n += 1

        self.prefix = text_output1.format(str(self.input_file_a), str(self.n), str(nukleotyd1), str(pozycja1+1))

    @selekcja
    def selekcja_sekwencji_2(self, pozycja1, nukleotyd1, pozycja2, nukleotyd2):
        """
        Funkcja filtrująca sekwencje po dwóch pozycjach nukleotydów.
        """
        for record in SeqIO.parse(self.input_seq, "fasta"):
            if record.seq[pozycja1].upper() == nukleotyd1 and\
               record.seq[pozycja2].upper() == nukleotyd2:
                self.output_seq += "{}\n{}\n".format(record.id, record.seq)
                self.n += 1

        self.prefix = text_output2.format(str(self.input_file_a), str(self.n),
                                          str(nukleotyd1), str(pozycja1+1),
                                          str(nukleotyd2), str(pozycja2+1))

    @selekcja
    def selekcja_sekwencji_3(self, pozycja1, nukleotyd1, pozycja2, nukleotyd2, pozycja3, nukleotyd3):
        """
        Funkcja filtrująca sekwencje po trzech pozycjach nukleotydów.
        """
        for record in SeqIO.parse(self.input_seq, "fasta"):
            if record.seq[pozycja1].upper() == nukleotyd1 and\
               record.seq[pozycja2].upper() == nukleotyd2 and\
               record.seq[pozycja3].upper() == nukleotyd3:
                self.output_seq += "{}\n{}\n".format(record.id, record.seq)
                self.n += 1

        self.prefix = text_output3 .format(str(self.input_file_a), str(self.n),
                                           str(nukleotyd1), str(pozycja1+1),
                                           str(nukleotyd2), str(pozycja2+1),
                                           str(nukleotyd3), str(pozycja3+1))

    def go(self):
        """
        Metoda wywołująca otwarcie pliku oraz filtrowanie sekencji. Metoda wywołuje również wiele innych metod.
        """
        self.open_input_file(self.input_file_a)
        self.input_control()
        self.get_pos_from_str()

        if self.nuc1 and self.nuc2 is False and self.nuc3 is False:
            self.selekcja_sekwencji_1(self.pos1, self.nuc1)

        elif self.nuc2 and self.nuc1 is False and self.nuc3 is False:
            self.selekcja_sekwencji_1(self.pos2, self.nuc2)

        elif self.nuc3 and self.nuc1 is False and self.nuc2 is False:
            self.selekcja_sekwencji_1(self.pos3, self.nuc3)

        elif self.nuc1 and self.nuc2 and self.nuc3 is False:
            self.selekcja_sekwencji_2(self.pos1, self.nuc1, self.pos2, self.nuc2)

        elif self.nuc1 and self.nuc3 and self.nuc2 is False:
            self.selekcja_sekwencji_2(self.pos1, self.nuc1, self.pos3, self.nuc3)

        elif self.nuc2 and self.nuc3 and self.nuc1 is False:
            self.selekcja_sekwencji_2(self.pos2, self.nuc2, self.pos3, self.nuc3)

        elif self.nuc1 and self.nuc2 and self.nuc3:
            self.selekcja_sekwencji_3(self.pos1, self.nuc1, self.pos2, self.nuc2, self.pos3, self.nuc3)

        else:
            self.output_seq = "brak sekwencji spełniających zadane warunki"

        self.text_field.delete('1.0', 'end')
        self.text_field.insert('1.0', self.output_seq)

    def interface_elem(self):
        """
        Tworzenie poszczególnych elementów GUI
        """

        szerokosc = 9

        label = Tkinter.Label(self.top, text="Aplikacja do filtrowania sekwencji", relief="ridge", pady=2, padx=300)
        label.grid(row=0, column=0, columnspan=szerokosc, sticky="nswe")

        # pusta linia oddzielająca okna input i out put
        Tkinter.Label(self.top).grid(row=1, column=0, columnspan=9)

        # szukanie pliku wejsciowego
        button = Tkinter.Button(self.top, text="Plik wejsciowy", command=self.input_file)
        button.grid(row=2, column=0, sticky="nswe")
        label = Tkinter.Label(self.top, width=7, textvariable=self.filepath_input)
        label.grid(row=2, column=1, columnspan=szerokosc-1, sticky="we")

        # pusta linia oddzielająca okna input i out put
        Tkinter.Label(self.top).grid(row=3, column=0, columnspan=szerokosc)

        # ComboBox 1
        label = Tkinter.Label(self.top, text=self.entry[0])
        label.grid(row=4, column=0)

        self.position1.grid(row=5, column=0, sticky="we")
        self.position1.insert(0, self.entry[0])
        self.position1['values'] = self.positions

        label = Tkinter.Label(self.top, text=self.entry[1])
        label.grid(row=4, column=1)

        self.nukleotyd1.grid(row=5, column=1, sticky="we")
        self.nukleotyd1.insert(0, self.entry[1])
        self.nukleotyd1['values'] = self.nukleotydy

        # ComboBox 2
        label = Tkinter.Label(self.top, text=self.entry[0])
        label.grid(row=4, column=3)

        self.position2.grid(row=5, column=3, sticky="we")
        self.position2.insert(0, self.entry[0])
        self.position2['values'] = self.positions

        label = Tkinter.Label(self.top, text=self.entry[1])
        label.grid(row=4, column=4)

        self.nukleotyd2.grid(row=5, column=4, sticky="we")
        self.nukleotyd2.insert(0, self.entry[1])
        self.nukleotyd2['values'] = self.nukleotydy

        # ComboBox 3
        label = Tkinter.Label(self.top, text=self.entry[0])
        label.grid(row=4, column=6)

        self.position3.grid(row=5, column=6, sticky="we")
        self.position3.insert(0, self.entry[0])
        self.position3['values'] = self.positions

        label = Tkinter.Label(self.top, text=self.entry[1])
        label.grid(row=4, column=7)

        self.nukleotyd3.grid(row=5, column=7, sticky="we")
        self.nukleotyd3.insert(0, self.entry[1])
        self.nukleotyd3['values'] = self.nukleotydy

        # pusta linia oddzielająca okna input i out put
        Tkinter.Label(self.top).grid(row=6, column=0, columnspan=szerokosc)

        button = Tkinter.Button(self.top, text="Let's Go !!", borderwidth=2,
                                command=self.go)
        button.grid(row=7, column=3, columnspan=2, sticky="nswe")

        # pusta linia oddzielająca okna input i out put
        Tkinter.Label(self.top).grid(row=8, column=0, columnspan=szerokosc)

        # pole tekstowe
        scrol = Tkinter.Scrollbar(self.top)
        scrol.grid(row=9, column=7, rowspan=3, sticky="ns")
        self.text_field = Tkinter.Text(self.top, width=123, yscrollcommand=scrol.set)
        self.text_field.grid(column=0, columnspan=8, row=9, rowspan=3)
        self.text_field.insert('1.0', "Tutaj będą wyświetlone wyniki przeszukiwania sekwencji")

        scrol.config(command=self.text_field.yview)

        # pusta linia
        Tkinter.Label(self.top).grid(row=10, column=0, columnspan=szerokosc)

        # Przycisk do zapisu danych
        button = Tkinter.Button(self.top, text="Save", borderwidth=2,
                                command=self.save)
        button.grid(row=14, column=3, columnspan=2, sticky="nswe")

        # pusta linia oddzielająca okna input i out put
        Tkinter.Label(self.top, text="Made by piemar@amu.edu.pl").grid(row=15, column=7)

        Tkinter.mainloop()

if __name__ == "__main__":
    prog = SequenceFinder()
    prog.interface_elem()
