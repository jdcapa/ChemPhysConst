"""A module which extracts Physical constants from an hdf5 database."""



import numpy as np
import sys
import os
import h5py

INT = np.int
FLOAT = np.float
here = os.path.dirname(__file__)
HDF5FILE = os.path.join(here, "data", "ChemPhysConst2016.hdf5")


class Constants(object):
    """
    Constants sourced from NIST.
    """
    def __init__(self):
        """Initiate the class, all info is provided through the HDF5 file."""
        super(Constants, self).__init__()
        self.angstroem = ["A", "Angs", "angs", "Angstroem", "angstroem"]
        self.constants = self.read_hdf5_data()

    def read_hdf5_data(self):
        """Read the physical_constants section of the provided HDF5 data."""
        if os.path.exists(HDF5FILE):
            constants = {}
            entries = ["description", "value", "uncertainty", "unit"]
            hdf5 = h5py.File(HDF5FILE, mode='r')
            for constant_name in list(hdf5["physical_constants"].keys()):
                const_attrs = hdf5["physical_constants"][constant_name].attrs
                constants[constant_name] = {e: const_attrs[e] for e in entries}
            hdf5.close()
            return constants
        else:
            sys.exit("The file {} is missing.".format(HDF5FILE))

    def nd(self, name):
        """Obtains the self.nist_data[name][value]"""
        return self.constants[name]['value']

    def planck_constant(self, unit="J*s"):
        """
        Return the Plack constant with a certain unit.

        Units:
            J*s
            eV*s
            Eh*s
        """
        if unit == "J*s":
            return self.nd('Planck_constant')
        elif unit == "eV*s":
            return self.nd('Planck_constant_in_eV_s')
        elif unit == "Eh*s":
            return self.nd('Planck_constant') * 1 / self.hartree_energy("J")

    def electron_mass(self, unit="u"):
        """
        Return the electron mass with a certain unit.

        Units:
            u
            kg
        """
        if unit == "kg":
            return self.nd('electron_mass')
        elif unit == "u":
            return self.nd('electron_mass_in_u')

    def neutron_mass(self, unit="u"):
        """Return the neutron mass with a certain unit.

        Units:
            u
            kg
        """
        if unit == "kg":
            return self.nd('neutron_mass')
        elif unit == "u":
            return self.nd('neutron_mass_in_u')

    def proton_mass(self, unit="u"):
        """
        Return the proton mass with a certain unit.

        Units:
            u
            kg
        """
        if unit == "kg":
            return self.nd('proton_mass')
        elif unit == "u":
            return self.nd('proton_mass_in_u')

    def avogadro_constant(self):
        """Return the Avogadro constant in mol^-1."""
        return self.nd('Avogadro_constant')

    def bohr_magneton(self):
        """Return the Bohr magneton in J/T."""
        return self.nd('Bohr_magneton')

    def bohr_radius(self, unit="m"):
        """Return the Bohr radius in m."""
        if unit == "m":
            return self.nd('Bohr_radius')
        elif unit in self.angstroem:
            return self.nd('Bohr_radius') * 1e10

    def Boltzmann_constant(self, unit="J/K"):
        """
        Return the Boltzmann constant with a certain unit.

        Units:
            J/K
            eV/K
            cm^-1/K
        """
        if unit == "J/K":
            return self.nd('Boltzmann_constant')
        elif unit == "eV/K":
            return self.nd('Boltzmann_constant_in_eV')
        elif unit == "cm^-1/K":
            k = self.nd('Boltzmann_constant_in_inverse_meters_per_kelvin')
            return k * 10e-2

    def hartree_energy(self, unit="kJ/mol"):
        """
        Return the Hartree energy with a certain unit.

        Units:
            J
            kJ/mol
            eV
        """
        if unit == "kJ/mol":
            Eh = self.nd('Hartree_energy')  # in J
            return Eh * self.avogadro_constant() * 1e-3
        elif unit == "J":
            return self.nd('Hartree_energy')
        elif unit == "kcal/mol":
            Eh = self.nd('Hartree_energy')  # in J
            return Eh * self.avogadro_constant() * 1e-3 / 4.184  # 4.184 J/cal
        elif unit == "eV":
            return self.nd('Hartree_energy_in_eV')

    def atomic_mass_constant(self, unit='kg'):
        """
        Return the atomic mass unit with a certain unit.

        Units:
            kg
            m_e (1 unit = x m_e)
        """
        if unit == 'kg':
            return self.nd('atomic_mass_constant')
        elif unit in ('me', 'm_e', 'atomic'):
            return 1.0 / self.electron_mass('u')

    def electron_volt(self):
        """Return eV in J."""
        return self.nd('electron_volt')

    def speed_of_light(self):
        """Return the speed of light in vacuum in m/s."""
        return self.nd('speed_of_light_in_vacuum')

    def standard_atmosphere(self):
        """Return the standard atmosphere in Pa."""
        return self.nd('standard_atmosphere')

    def standard_temperature(self):
        """Return the standard temperature in K."""
        return 298.15

    def molar_gas_constant(self):
        """Return the molar gas constant R in J mol^-1 K^-1."""
        return self.nd('molar_gas_constant')

    def third_order_LeviCevita_Tensor(self):
        """
        Return a third order Levi_Civita tensor e_ijk.

        Here,
         eijk[0, 1, 2] = eijk[1, 2, 0] = eijk[2, 0, 1] = 1
         eijk[0, 2, 1] = eijk[2, 1, 0] = eijk[1, 0, 2] = -1
        """
        e = np.array(
            [
                [
                    [int((i - j) * (j - k) * (k - i) / 2) for k in range(3)]
                    for j in range(3)]
                for i in range(3)
            ], dtype=FLOAT)
        return e
