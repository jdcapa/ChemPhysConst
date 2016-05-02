# -*- coding: UTF-8 -*-
"""This module contains all the nuclide and element specific classes."""

from __future__ import print_function
from __future__ import division
# from __future__ import unicode_literals  # This breaks the rb+ writing (hdf5)
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
DATAFILE = os.path.join(here, "data", "2015_mass_weight_abundance.dat")
ELEMFILE = os.path.join(here, "data", "elements.csv")
HDF5FILE = os.path.join(here, "data", "2015_mass_weight_abundance.hdf5")


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

    Variables: symbol, number, representing_isotope, name, type, group, period
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


class PeriodicTableNIST(object):
    """
    Reads the Periodic Table information from the provided NIST DATAFILE.

    Furthermore it gets additional element data from a self-curated
     'elements.csv' file.
    This is slow and ugly! Hence this is only used for the initial database
    creation.
      Variables: an_to_sy, sy_to_an
      Methods: element(mass_number or symbol)
    """

    def __init__(self):
        """Initiate the class, all info is provided through files."""
        super(PeriodicTableNIST, self).__init__()
        self.an_to_sy = self.atomic_number_to_symbol()
        self.sy_to_an = self.symbol_to_atomic_number()
        self.element_properties = self.interpret_properties()

    def element(self, element):
        """
        Return the corresponding class to an element.

        The variable element can be int (atomic number) or str (symbol).
        """
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

        isotopes, aw_str = self.element_isotopic_data(an)
        kwargs = self.element_properties[an]
        kwargs["isotopes"] = isotopes
        kwargs["atomic_weight_str"] = aw_str
        return Element(sy, an, **kwargs)

    def atomic_number_to_symbol(self):
        """Return a dictionary mapping the atomic numbers to the symbols."""
        re_atomic_number = re.compile('Atomic Number = (\d+)')
        re_atomic_symbol = re.compile('Atomic Symbol = (\w+)')
        element_map = []
        check = []
        for line in self.read_data_file():
            if re_atomic_number.search(line):
                num = INT(re_atomic_number.search(line).group(1))
            if re_atomic_symbol.search(line):
                sy = re_atomic_symbol.search(line).group(1)
                if num not in check:
                    check.append(num)
                    element_map.append([num, sy])
        element_map.sort(key=operator.itemgetter(0))
        return OrderedDict(element_map)

    def symbol_to_atomic_number(self):
        """Return a dictionary mapping the the atomic symbols to numbers."""
        return OrderedDict((sy, an) for an, sy in self.an_to_sy.items())

    def element_isotopic_data(self, atomic_num):
        """
        Return a list containing all the isotopes for an atomic_number.

        This in turn is Nuclide class containing masses and abundances.
        Also returns the NIST standard weight string.
        """
        re_atomic_number = re.compile('Atomic Number = (\d+)')
        re_atomic_symbol = re.compile('Atomic Symbol = (\w+)')
        re_mass_number = re.compile('Mass Number = (\d+)')
        re_atomic_mass = re.compile('Relative Atomic Mass = ([\d.]+)')
        re_isotopic_abundance = re.compile('Isotopic Composition = ([\d.]*)')
        re_atomic_weight = re.compile('Standard Atomic Weight = ([\[\d.,\]]+)')
        isotopes = []
        mass_number = 0
        atomic_mass = 0.0
        atomic_symbol = ''
        isotopic_abund = 0.0
        atomic_weight_str = ''
        read_flag = False
        for line in self.read_data_file():
            # Searching for the right element
            if re_atomic_number.search(line):
                if INT(re_atomic_number.search(line).group(1)) == atomic_num:
                    if mass_number:
                        isotope = Nuclide(atomic_symbol, atomic_num,
                                          mass_number,
                                          atomic_mass=atomic_mass,
                                          abundance=isotopic_abund)
                        isotopes.append(isotope)
                    read_flag = True
                    mass_number = 0
                    atomic_mass = 0.0
                    atomic_symbol = ''
                    isotopic_abund = 0.0
                    continue
            if not read_flag:
                continue
            # Searching for properties
            elif re_atomic_symbol.search(line):
                atomic_symbol = re_atomic_symbol.search(line).group(1)
                continue
            elif re_mass_number.search(line):
                mass_number = INT(re_mass_number.search(line).group(1))
                continue
            elif re_atomic_mass.search(line):
                atomic_mass = FLOAT(re_atomic_mass.search(line).group(1))
                continue
            elif re_isotopic_abundance.search(line):
                tmp_isoab = re_isotopic_abundance.search(line).group(1)
                if tmp_isoab:
                    isotopic_abund = FLOAT(tmp_isoab)
                else:
                    isotopic_abund = 0.0
                continue
            elif (re_atomic_weight.search(line)):
                if not atomic_weight_str:
                    atomic_weight_str = re_atomic_weight.search(line).group(1)
                read_flag = False
                continue

        return isotopes, atomic_weight_str

    def read_data_file(self):
        """
        Return the data file as an string.

        The file is always small, so this is ok and makes everything faster.
        """
        with open(DATAFILE) as data_raw:
            data = data_raw.read()
        return data.split('\n')

    def read_properties_file(self):
        """Return the comma separated data from the ELEMFILE."""
        with open(ELEMFILE) as data_raw:
            data = data_raw.read()
        return [line.strip().split(',') for line in data.split('\n')[2:]
                if line.strip().split(',')]

    def interpret_properties(self):
        """Interpret the property list and return an OrderedDict."""
        def check_type(entry, convert_type):
            try:
                return convert_type(entry)
            except:
                return None

        raw_prop = self.read_properties_file()
        prop = OrderedDict()
        for line in raw_prop:
            if len(line) < 17:
                continue
            z = check_type(line[0], int)
            el_type = check_type(line[1], str)
            name = check_type(line[3], str)
            group = check_type(line[4], int)
            periode = check_type(line[5], int)
            density = check_type(line[6], float)
            melting_point = check_type(line[7], float)
            boiling_point = check_type(line[8], float)
            electro_negativity = check_type(line[10], float)
            abundance_crust = check_type(line[11], float)
            covalent_r_single = check_type(line[13], int)
            covalent_r_double = check_type(line[14], int)
            covalent_r_triple = check_type(line[15], int)
            vdW_r = check_type(line[16], int)
            e_config = check_type(line[17], str)

            propset = {"name": name,
                       "element_type": el_type,
                       "group": group,
                       "periode": periode,
                       "density": density,
                       "melting_point": melting_point,
                       "boiling_point": boiling_point,
                       "electro_negativity": electro_negativity,
                       "abundance_crust": abundance_crust,
                       "covalent_r_single": covalent_r_single,
                       "covalent_r_double": covalent_r_double,
                       "covalent_r_triple": covalent_r_triple,
                       "vdW_r": vdW_r,
                       "e_config": e_config}
            prop[z] = {k: v for k, v in propset.items() if v}
        return prop


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

    def __del__(self):
        """Make sure the HDF5 is closed upon cleanup."""
        self.hdf5.close()

    def atomic_number_to_symbol(self):
        """Return a dictionary mapping the atomic numbers to the symbols."""
        el_map = [[el.attrs['number'], sy] for sy, el in self.hdf5.items()]
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
                      "abundance_crust", "covalent_r_single",
                      "covalent_r_double", "covalent_r_triple", "vdW_r",
                      "e_config"]
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

        el = self.hdf5[sy].attrs
        el_cp = {k: el[k] for k in properties}
        isotopes = []
        # print(self.hdf5[sy]['2'].attrs['atomic_symbol'])
        # sys.exit()
        for mn in self.hdf5[sy].keys():
            iso = self.hdf5[sy][mn].attrs
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


class ExportPeriodicTableNIST(object):
    """This class harbours the export functions for the PeriodicTable class."""

    def __init__(self):
        """Initiate the class, no input, just export."""
        super(ExportPeriodicTableNIST, self).__init__()
        self.PTN = PeriodicTableNIST()

    def print_masses_cpp_array(self):
        """
        Print the representative masses for export to an ugly C++ array.

        I'm looking at you Orca!
        """
        header = "    //{0:^10}       {1:^3}  {2:<4}   {3:^12}\n".format(
            "m[a]", "Z", "Sy", "<m[a]>")
        for an, sy in self.PTN.an_to_sy.items():
            elem = self.PTN.element(an)
            tmp_str = "    {:<15}".format(str(elem.mass) + ",")
            tmp_str += " // {:03d}".format(an)
            tmp_str += "  {:<4}".format(sy)
            tmp_str += " {:>12,.8f}\n".format(elem.weight)
            header += tmp_str
        print(header)

    def export_elements_to_hdf5(self):
        """
        Write a hdf5 file dumping the elements data.

        This should speed up subsequent loading of the element data.
        We are dumping:
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
        userblock_size = 512
        user_block_string = ("PeriodicTable data created from the \n" +
                             "Atomic Weights and Isotopic Compositions for " +
                             "All Elements compilation (2015)\n" +
                             "http://www.nist.gov/pml/data/comp.cfm\n" +
                             "Created by the PeriodicTable Library:\n" +
                             "https://github.com/jdcapa/PeriodicTable\n")

        if userblock_size < len(user_block_string):
            sys.exit('userblock_size is too small')
        hdf5 = h5py.File(HDF5FILE, mode='w', userblock_size=userblock_size)
        hdf5 = self.write_header(hdf5)

        for an, sy in self.PTN.an_to_sy.items():
            el = self.PTN.element(an)
            el_grp = hdf5.create_group(sy)
            el_grp.attrs.create('symbol', sy)
            el_grp.attrs.create('number', an)
            el_grp.attrs.create('name', el.name)
            el_grp.attrs.create('atomic_weight_str', el.atomic_weight_str)
            el_grp.attrs.create('element_type', el.element_type)
            el_grp.attrs.create('atomic_weight_str', el.atomic_weight_str)
            el_grp.attrs.create('density', el.density)
            el_grp.attrs.create('melting_point', el.melting_point)
            el_grp.attrs.create('boiling_point', el.boiling_point)
            el_grp.attrs.create('electro_negativity', el.electro_negativity)
            el_grp.attrs.create('abundance_crust', el.abundance_crust)
            el_grp.attrs.create('covalent_r_single', el.covalent_r_single)
            el_grp.attrs.create('covalent_r_double', el.covalent_r_double)
            el_grp.attrs.create('covalent_r_triple', el.covalent_r_triple)
            el_grp.attrs.create('vdW_r', el.vdW_r)
            el_grp.attrs.create('e_config', el.e_config)
            for isotope in el.isotopes:
                iso_grp = el_grp.create_group(str(isotope.mass_number))
                iso_grp.attrs.create('atomic_symbol', isotope.atomic_symbol)
                iso_grp.attrs.create('atomic_number', isotope.atomic_number)
                iso_grp.attrs.create('mass_number', isotope.mass_number)
                iso_grp.attrs.create('atomic_mass', isotope.atomic_mass)
                iso_grp.attrs.create('abundance', isotope.abundance)

        hdf5.flush()
        hdf5.close()
        with open(HDF5FILE, 'rb+') as hdf5_ub:
            hdf5_ub.write(user_block_string)

    def write_header(self, f):
        """Add some root attributes to an hdf5 object."""
        from time import gmtime, strftime

        f.attrs['file_name'] = os.path.basename(HDF5FILE)
        f.attrs['file_time'] = strftime("%Y-%m-%d %H:%M:%S", gmtime())
        f.attrs['creator'] = 'periodictable.py'
        f.attrs['HDF5_Version'] = h5py.version.hdf5_version
        f.attrs['h5py_version'] = h5py.version.version
        return f
