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

    def read_data_file(self):
        """
        Returns the data file as an string.
        The file is always small, so this is ok and makes everything faster.
        """
        with open(DATAFILE) as data_raw:
            data = data_raw.read()
        return data.split('\n')
