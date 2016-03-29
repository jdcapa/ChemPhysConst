#!/usr/bin/env python


import re
import os
# from __future__ import division
# from __future__ import print_function
# from __future__ import unicode_literals

DATAFILE = os.path.join("data", "2015_mass_weight_abundance.dat")
# DATAFILE = "2015_mass_weight_abundance.dat"


class PeriodicTable(object):
    """PeriodicTable"""
    def __init__(self):
        super(PeriodicTable, self).__init__()
        self.data = self.read_data_file()
        self.atomic_numbers = self.element_map()

    def element_map(self):
        """
        Returns a dictionary mapping the atomic symbols to the atomic
        """
        re_atomic_number = re.compile('Atomic Number = (\d+)')
        re_atomic_symbol = re.compile('Atomic Symbol = (\w+)')
        element_map = {}
        for line in self.data:
            if re_atomic_number.search(line):
                num = int(re_atomic_number.search(line).group(1))
            if re_atomic_symbol.search(line):
                sym = re_atomic_symbol.search(line).group(1)
                if num not in element_map:
                    element_map[num] = sym
        return element_map

    def element_isotopic_data(self, atomic_num):
        '''
        Returns a dictionary for atomic_number containing all the isotopes
         which in turn are dictionaries containing masses and abundances.
        '''
        re_atomic_number = re.compile('Atomic Number = (\d+)')
        re_mass_number = re.compile('Mass Number = (\d+)')
        re_atomic_mass = re.compile('Relative Atomic Mass = ([\d.]+)')
        re_isotopic_abundance = re.compile('Isotopic Composition = ([\d.]*)')
        re_atomic_weight = re.compile('Standard Atomic Weight = ([\[\d.,\]]+)')
        isotopes = {}
        read_flag = False
        standard_status = -1
        for line in self.data:
            # Searching for the right element
            if re_atomic_number.search(line):
                if int(re_atomic_number.search(line).group(1)) == atomic_num:
                    read_flag = True
                    mass_number = 0
                    atomic_mass = 0.0
                    isotopic_abundance = 0.0
                    atomic_weight = 0.0
                    continue
            if not read_flag:
                continue

            # Searching for properties
            if re_mass_number.search(line):
                mass_number = int(re_mass_number.search(line).group(1))
                continue
            if re_atomic_mass.search(line):
                atomic_mass = float(re_atomic_mass.search(line).group(1))
                continue
            if re_isotopic_abundance.search(line):
                ab_str = re_isotopic_abundance.search(line).group(1)
                if ab_str:
                    isotopic_abundance = float(ab_str)
                else:
                    isotopic_abundance = 0.0
                continue

            # Element properties
            if standard_status > -1:
                if (re_mass_number.search(line)):
                    pass
                if (re_mass_number.search(line)):
                    pass

    def read_data_file(self):
        """
        Returns the data file as an string.
        The file is always small, so this is ok and makes everything faster.
        """
        with open(DATAFILE) as data_raw:
            data = data_raw.read()
        return data.split('\n')
