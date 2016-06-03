# -*- coding: UTF-8 -*-
"""This module contains all the nuclide and element specific classes."""

import sys
import re
import os
from collections import OrderedDict
import operator
import h5py
import numpy as np

INT = np.int
FLOAT = np.float
here = os.path.dirname(__file__)
HDF5FILE = os.path.join(here, "data", "ChemPhysConst2016.hdf5")


class Nuclide(object):
    """
    This contains properties of a given nuclide.

    Var: atomic_symbol, atomic_number, mass_number, atomic_mass, abundance
    """

    def __init__(self, sy, an, mn, **kwargs):
        """Initiate symbol, atomic number, mass number and kewords."""
        super(Nuclide, self).__init__()
        self.atomic_symbol = sy
        self.atomic_number = an
        self.mass_number = mn
        self.nuclide_setter(kwargs)
        self.stability()

    def nuclide_setter(self, kwargs):
        """Interpret kwargs (adds atomic mass, isotopic abundance)."""
        if "atomic_mass" in kwargs:
            self.atomic_mass = kwargs["atomic_mass"]
        else:
            self.atomic_mass = FLOAT(0)
        if "abundance" in kwargs:
            self.abundance = kwargs["abundance"]
        else:
            self.abundance = FLOAT(0)

    def stability(self):
        """Check if the isotope is stable."""
        self.is_stable = False
        if self.abundance > 0:
            self.is_stable = True


class Element(object):
    """
    This contains properties of a given element.

    Variables: symbol, number, representing_isotope, name, type
    Properties: is_stable (True|False)
                mass (in u)
                weight (in u)
                density (in g/cm^3)
                melting_point (in K)
                boiling_point (in K)
                electro_negativity (Pauling scale)
                abundance_crust (in mg/kg)
                covalent_r_single (in pm)
                covalent_r_double (in pm)
                covalent_r_triple (in pm)
                vdW_r (in pm)
                e_config (html string)
                group
                period

    Methods: isotope(mass_number)
    """

    def __init__(self, sy, an, **kwargs):
        """Initiate symbol, atomic number and kewords."""
        super(Element, self).__init__()
        self.symbol = sy
        self.number = an
        self.element_setter(kwargs)
        self.is_stable = self.stability()
        self.representing_isotope = self.get_representative_isotope()
        self.mass = self.representing_isotope.atomic_mass
        self.weight = self.calculate_atomic_weight()

    def element_setter(self, kwargs):
        """
        Interpret kwargs.

        Adds atomic mass, isotopic abundance and if it is stable.
        """
        def check_kwargs_for(keyword, default):
            if keyword in kwargs:
                if kwargs[keyword]:
                    return kwargs[keyword]
            return default
        self.isotopes = check_kwargs_for("isotopes", [])
        self.atomic_weight_str = check_kwargs_for("atomic_weight_str", "")
        self.name = check_kwargs_for("name", "")
        self.element_type = check_kwargs_for("element_type", "unknown")
        self.group = check_kwargs_for("group", 0)
        self.period = check_kwargs_for("period", 0)
        self.density = check_kwargs_for("density", 0.0)
        self.melting_point = check_kwargs_for("melting_point", 0.0)
        self.boiling_point = check_kwargs_for("boiling_point", 0.0)
        self.electro_negativity = check_kwargs_for("electro_negativity", 0.0)
        self.abundance_crust = check_kwargs_for("abundance_crust", 0.0)
        # Covalent Bond distances according to Pyykkö
        self.covalent_r_single = check_kwargs_for("covalent_r_single", 150.0)
        self.covalent_r_double = check_kwargs_for("covalent_r_double", 0.0)
        self.covalent_r_triple = check_kwargs_for("covalent_r_triple", 0.0)
        # Misc
        self.vdW_r = check_kwargs_for("vdW_r", 2.00)
        self.e_config = check_kwargs_for("e_config", "")

    def isotope(self, mass_number):
        """Return an nuclide object of a given isotope."""
        for isotope in self.isotopes:
            if isotope.mass_number == mass_number:
                return isotope
        sys.exit("Isotope {}-{} not known.".format(self.atomic_symbol,
                                                   mass_number))

    def stability(self):
        """Check if the element is stable."""
        for isotope in self.isotopes:
            if isotope.abundance > 0.0:
                return True
        return False

    def calculate_atomic_weight(self, mode='average'):
        """
        Return the atomic weight.

        This is either a weighted average over the stable isotopes
         (mode='average'), the average of the range given in the data file
         (mode='range') or in case of an unstable element the mass of the most
         stable isotope.
        """
        if self.is_stable:
            if mode == 'average':
                atomic_weight = 0.0
                for isotope in self.isotopes:
                    atomic_weight += isotope.abundance * isotope.atomic_mass
                if atomic_weight:
                    return atomic_weight
                else:
                    return self.calculate_atomic_weight('range')
            elif mode == 'range':
                if ',' in self.atomic_weight_str:
                    # If there is a range, we should assign the average
                    re_num = re.compile('\[(\d.+),(\d.+)\]')
                    masses = [re_num.match(self.atomic_weight_str).group(1),
                              re_num.match(self.atomic_weight_str).group(2)]
                    atomic_weight = (FLOAT(masses[0]) + FLOAT(masses[1])) / 2
                else:
                    try:
                        # This is just a number
                        atomic_weight = FLOAT(self.atomic_weight_str)
                    except:
                        sys.exit("Couldn't convert the standard atomic Weight")
                return atomic_weight
        else:
            return self.representing_isotope.atomic_mass

    def get_representative_isotope(self):
        """
        Return a representing isotope (nuclide class).

        Depending on whether or not the element is stable, the representing
         isotope is either the most abundant or the one with the longest
         half-life.
        """
        if self.is_stable:
            check = [0, 0.0]
            for isotope in self.isotopes:
                if isotope.abundance > check[1]:
                    check[0] = isotope.mass_number
                    check[1] = isotope.abundance
            return self.isotope(check[0])
        else:
            ma = INT(re.match('\[(\d+)\]', self.atomic_weight_str).group(1))
            masses = dict((isotope.mass_number, isotope.atomic_mass)
                          for isotope in self.isotopes)
            diff = [0, 100.00]
            for mn, mass in masses.items():
                if abs(mass - ma) < diff[1]:
                    diff[0] = mn
                    diff[1] = abs(mass - ma)
            return self.isotope(diff[0])


class PeriodicTable(object):
    """Reads the Periodic Table information from the provided HDF5 DATAFILE."""

    def __init__(self):
        """Initiate the class, all info is provided through the HDF5 file."""
        super(PeriodicTable, self).__init__()
        if os.path.exists(HDF5FILE):
            self.hdf5 = h5py.File(HDF5FILE, mode='r')
            self.an_to_sy = self.atomic_number_to_symbol()
            self.sy_to_an = self.symbol_to_atomic_number()
        else:
            sys.exit("The file {} is missing.".format(HDF5FILE))

    # def __del__(self):
    #     """Make sure the HDF5 is closed upon cleanup."""
    #     self.hdf5.close()

    def atomic_number_to_symbol(self):
        """Return a dictionary mapping the atomic numbers to the symbols."""
        el_map = [[el.attrs['number'], sy] for sy, el in
                  self.hdf5["periodic_table"].items()]
        el_map.sort(key=operator.itemgetter(0))
        return OrderedDict(el_map)

    def symbol_to_atomic_number(self):
        """Return a dictionary mapping the the atomic symbols to numbers."""
        return OrderedDict((sy, an) for an, sy in self.an_to_sy.items())

    def element(self, element):
        """
        Import the generated hdf5 version of the Periodic table data.

        We are reading:
            Nuclides:
                atomic_symbol:              str
                atomic_number:              int
                mass_number:                int
                atomic_mass:                float
                abundance:                  float
            Elements:
                symbol:                     str
                number:                     int
                name:                       str
                element_type:               str
                atomic_weight_str:          str
                isotopes:                   Nuclides()
                group                       int
                period                      int
                density:                    float       (in g/cm^3)
                melting_point:              float       (in K)
                boiling_point:              float       (in K)
                electro_negativity:         float       (Pauling scale)
                abundance_crust:            float       (in mg/kg)
                covalent_r_single:          int         (in pm)
                covalent_r_double:          int         (in pm)
                covalent_r_triple:          int         (in pm)
                vdW_r:                      int         (in pm)
                e_config:                   str         (html string)
        """
        properties = ["name", "element_type", "atomic_weight_str", "density",
                      "melting_point", "boiling_point", "electro_negativity",
                      "group", "period", "abundance_crust",
                      "covalent_r_single", "covalent_r_double",
                      "covalent_r_triple", "vdW_r", "e_config"]
        try:
            if type(element) is INT:
                sy = self.an_to_sy[element]
                an = element
            elif type(element) is str:
                sy = element
                an = self.sy_to_an[element]
        except KeyError as e:
            raise e("get_element_properties() requires either int or str type."
                    " {} was provided.".format(type(element)))

        el = self.hdf5["periodic_table"][sy].attrs
        el_cp = {k: el[k] for k in properties}
        isotopes = []
        # print(self.hdf5[sy]['2'].attrs['atomic_symbol'])
        # sys.exit()
        for mn in self.hdf5["periodic_table"][sy].keys():
            iso = self.hdf5["periodic_table"][sy][mn].attrs
            mass_number = INT(mn)
            atomic_symbol = iso['atomic_symbol']
            atomic_num = iso['atomic_number']
            atomic_mass = iso['atomic_mass']
            isotopic_abund = iso['abundance']
            isotope = Nuclide(atomic_symbol, atomic_num, mass_number,
                              atomic_mass=atomic_mass,
                              abundance=isotopic_abund)
            isotopes.append(isotope)
        el_cp["isotopes"] = isotopes
        return Element(sy, an, **el_cp)
