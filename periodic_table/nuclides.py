#!/usr/bin/env python


import re
import os
# from __future__ import division
# from __future__ import print_function
# from __future__ import unicode_literals

DATAFILE = os.path.join("data", "2015_mass_weight_abundance.dat")


class PeriodicTable(object):
    """PeriodicTable"""
    def __init__(self):
        super(PeriodicTable, self).__init__()
        self.data = self.read_data_file()

    def atomic_number_to_symbol(self):
        """
        Returns a dictionary mapping the atomic numbers to  the symbols.
        """
        re_atomic_number = re.compile('Atomic Number = (\d+)')
        re_atomic_symbol = re.compile('Atomic Symbol = (\w+)')
        element_map = {}
        for line in self.data:
            if re_atomic_number.search(line):
                num = int(re_atomic_number.search(line).group(1))
            if re_atomic_symbol.search(line):
                sy = re_atomic_symbol.search(line).group(1)
                if num not in element_map:
                    element_map[num] = sy
        return element_map

    def symbol_to_atomic_number(self):
        """
        Returns a dictionary mapping the the atomic symbols to the numbers.
        """
        an_to_sy = self.atomic_number_to_symbol()
        return dict((sy, an) for an, sy in an_to_sy.items())

    def element_isotopic_data(self, atomic_num):
        '''
        Returns a dictionary for atomic_number containing all the isotopes
         which in turn are dictionaries containing masses and abundances.
        '''
        re_atomic_number = re.compile('Atomic Number = (\d+)')
        re_atomic_symbol = re.compile('Atomic Symbol = (\w+)')
        re_mass_number = re.compile('Mass Number = (\d+)')
        re_atomic_mass = re.compile('Relative Atomic Mass = ([\d.]+)')
        re_isotopic_abundance = re.compile('Isotopic Composition = ([\d.]*)')
        re_atomic_weight = re.compile('Standard Atomic Weight = ([\[\d.,\]]+)')
        isotopes = {}
        mass_number = 0
        atomic_mass = 0.0
        atomic_symbol = ''
        isotopic_abundance = 0.0
        tmp_aw = ''
        read_flag = False
        standard_status = -1
        for line in self.data:
            # Searching for the right element
            if re_atomic_number.search(line):
                if int(re_atomic_number.search(line).group(1)) == atomic_num:
                    if mass_number:
                        isotopes[mass_number] = {
                            'atomic number': atomic_num,
                            'atomic symbol': atomic_symbol,
                            'atomic mass': atomic_mass,
                            'isotopic abundance': isotopic_abundance}
                    read_flag = True
                    mass_number = 0
                    atomic_mass = 0.0
                    atomic_symbol = ''
                    isotopic_abundance = 0.0
                    continue
            if not read_flag:
                continue
            # Searching for properties
            if re_atomic_symbol.search(line):
                atomic_symbol = re_atomic_symbol.search(line).group(1)
                continue
            if re_mass_number.search(line):
                mass_number = int(re_mass_number.search(line).group(1))
                continue
            if re_atomic_mass.search(line):
                atomic_mass = float(re_atomic_mass.search(line).group(1))
                continue
            if re_isotopic_abundance.search(line):
                tmp_isoab = re_isotopic_abundance.search(line).group(1)
                if tmp_isoab:
                    isotopic_abundance = float(tmp_isoab)
                else:
                    isotopic_abundance = 0.0
                continue
            # Element properties
            if (re_atomic_weight.search(line)):
                if standard_status == -1:
                    tmp_aw = re_atomic_weight.search(line).group(1)
                    standard_status = 0
                read_flag = False
                continue
        if tmp_aw:
            isotopes = self.atomic_weight(isotopes, tmp_aw)
        return isotopes

    def atomic_weight(self, isotopes, tmp_aw):
        '''
        Adds an 'atomic weight' and 'stable' key to the isotopes dict.
        '''
        if ',' in tmp_aw:
            # If there is a range, we should assign the average
            re_num = re.compile('\[(\d.+),(\d.+)\]')
            masses = [re_num.match(tmp_aw).group(1),
                      re_num.match(tmp_aw).group(2)]
            isotopes['atomic weight'] = (float(masses[0]) +
                                         float(masses[1])) / 2
            isotopes['stable'] = True
        elif '[' in tmp_aw:
            # We should assign the mass of the most stable isotope
            masses = [value['atomic mass'] for value in isotopes.values()]
            ma = int(re.match('\[(\d+)\]', tmp_aw).group(1))
            isotopes['atomic weight'] = min(masses, key=lambda x: abs(x - ma))
            isotopes['stable'] = False
        else:
            try:
                # This is just a number
                isotopes['atomic weight'] = float(tmp_aw)
                isotopes['stable'] = True
            except:
                print 'Could not convert the Standard Atomic Weight'
        return isotopes

    def read_data_file(self):
        """
        Returns the data file as an string.
        The file is always small, so this is ok and makes everything faster.
        """
        with open(DATAFILE) as data_raw:
            data = data_raw.read()
        return data.split('\n')
