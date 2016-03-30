#!/usr/bin/env python


import re
import os
import sys
# from __future__ import division
# from __future__ import print_function
# from __future__ import unicode_literals

DATAFILE = os.path.join("data", "2015_mass_weight_abundance.dat")


class PeriodicTable(object):
    """PeriodicTable"""

    def __init__(self):
        super(PeriodicTable, self).__init__()
        self.data = self.read_data_file()

    def get_elements(self):
        """
        Constructs the element dictionary
        """
        elements = {}
        sy_to_an = self.symbol_to_atomic_number()
        for sy, an in sy_to_an.items():
            elements[sy] = {}
            isotopes = self.element_isotopic_data(an)
            atomic_weight = self.atomic_weight(isotopes, 'average')
            element_is_stable = self.elemtent_stability(isotopes)
            representing_isotope = self.get_representative_isotope(isotopes)
            representing_mass = isotopes[representing_isotope]['atomic mass']

            elements[sy]['atomic number'] = an
            elements[sy]['atomic symbol'] = sy
            elements[sy]['atomic weight'] = atomic_weight
            elements[sy]['is stable'] = element_is_stable
            elements[sy]['isotopes'] = isotopes
            elements[sy]['representing isotope'] = representing_isotope
            elements[sy]['representing mass'] = representing_mass
        return elements

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
        atomic_weight = ''
        read_flag = False
        for line in self.data:
            # Searching for the right element
            if re_atomic_number.search(line):
                if int(re_atomic_number.search(line).group(1)) == atomic_num:
                    if mass_number:
                        isotopes[mass_number] = {
                            'atomic number': atomic_num,
                            'atomic symbol': atomic_symbol,
                            'atomic mass': atomic_mass,
                            'isotopic abundance': isotopic_abundance,
                            'atomic weight': atomic_weight}
                    read_flag = True
                    mass_number = 0
                    atomic_mass = 0.0
                    atomic_symbol = ''
                    isotopic_abundance = 0.0
                    continue
            if not read_flag:
                continue
            # Searching for properties
            elif re_atomic_symbol.search(line):
                atomic_symbol = re_atomic_symbol.search(line).group(1)
                continue
            elif re_mass_number.search(line):
                mass_number = int(re_mass_number.search(line).group(1))
                continue
            elif re_atomic_mass.search(line):
                atomic_mass = float(re_atomic_mass.search(line).group(1))
                continue
            elif re_isotopic_abundance.search(line):
                tmp_isoab = re_isotopic_abundance.search(line).group(1)
                if tmp_isoab:
                    isotopic_abundance = float(tmp_isoab)
                else:
                    isotopic_abundance = 0.0
                continue
            elif (re_atomic_weight.search(line)):
                atomic_weight = re_atomic_weight.search(line).group(1)
                read_flag = False
                continue

        return isotopes

    def elemtent_stability(self, isotopes):
        '''
        Finds out whether or not an element is stable.
        '''
        aw_str = isotopes[isotopes.keys()[0]]['atomic weight']
        if not aw_str:
            return False
        elif '[' in aw_str:
            if ',' not in aw_str:
                return False
            else:
                return True
        else:
            return True

    def atomic_weight(self, isotopes, mode='average'):
        """
        Returns the atomic weight which is either a weighted average over the
         stable isotopes (mode='average'), the average of the range given in
         the data file (mode='range') or in case of an unstable element the
         mass of the most stable isotope.
        """
        element_is_stable = self.elemtent_stability(isotopes)
        if element_is_stable:
            if mode == 'average':
                atomic_weight = 0.0
                for properties in isotopes.values():
                    atomic_weight += (properties['isotopic abundance'] *
                                      properties['atomic mass'])
                if atomic_weight:
                    return atomic_weight
                else:
                    return self.atomic_weight(isotopes, 'range')
            elif mode == 'range':
                tmp_aw = isotopes[isotopes.keys()[0]]['atomic weight']
                if ',' in tmp_aw:
                    # If there is a range, we should assign the average
                    re_num = re.compile('\[(\d.+),(\d.+)\]')
                    masses = [re_num.match(tmp_aw).group(1),
                              re_num.match(tmp_aw).group(2)]
                    atomic_weight = (float(masses[0]) + float(masses[1])) / 2
                else:
                    try:
                        # This is just a number
                        atomic_weight = float(tmp_aw)
                    except:
                        sys.exit("Couldn't convert the standard atomic Weight")
                return atomic_weight
        else:
            rep_isotope = self.get_representative_isotope(isotopes)
            return isotopes[rep_isotope]['atomic mass']

    def get_representative_isotope(self, isotopes):
        """
        Depending on whether or not the element is stable, the representing
         isotope is either the most abundant or the one with the longest
         half-life.
        Returns a mass number.
        """
        element_is_stable = self.elemtent_stability(isotopes)
        if element_is_stable:
            check = [0, 0.0]
            for mass_number, properties in isotopes.items():
                if properties['isotopic abundance'] > check[1]:
                    check[0] = mass_number
                    check[1] = properties['isotopic abundance']
            return check[0]
        else:
            tmp_aw = isotopes[isotopes.keys()[0]]['atomic weight']
            ma = int(re.match('\[(\d+)\]', tmp_aw).group(1))
            masses = dict((mn, v['atomic mass']) for mn, v in isotopes.items())
            diff = [0, 100.00]
            for mn, mass in masses.items():
                if abs(mass - ma) < diff[1]:
                    diff[0] = mn
                    diff[1] = abs(mass - ma)
            return diff[0]

    def read_data_file(self):
        """
        Returns the data file as an string.
        The file is always small, so this is ok and makes everything faster.
        """
        with open(DATAFILE) as data_raw:
            data = data_raw.read()
        return data.split('\n')
