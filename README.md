# ChemPhysConst


A lightweight python library allowing access to physical constants, periodic table
 data and mathematical conventions.
The periodic table data contains atomic weights, isotopic compositions and
 elemental properties for all known chemical elements.

## Installation

**Requirements:** Python 2.7, pip,
 [hd5py](http://docs.h5py.org/en/latest/index.html)

Simply use pip to install the periodic table python library:

    $ pip install --user git+https://github.com/jdcapa/ChemPhysConst

## Usage

The following briefly describes the import and scope of the `PeriodicTable`
 and `Constants` module.

### PeriodicTable
-----

In your python import section put:

```python
from chemphysconst import PeriodicTable
```

Initialise:

```python
PT = PeriodicTable()
```

where PT contains the dictionaries

```python
PT.an_to_sy               # mapping all the atomic numbers to the symbols
PT.sy_to_an               # mapping all the symbols to the atomic numbers
```

Access elements via the PT.elements() method:

```python
el_H = PT.element(1)      # Either via the atomic number (Hydrogen)
el_F = PT.element('F')    # or the atomic symbol (Fluorine)
el_Tc = PT.element(43)    # Technetium
el_U = PT.element('U')    # Uranium
```

The elements contains loads of properties themselves:

**Basic properties**:

```python
el_H.symbol               # atomic number
el_F.number               # atomic symbol
el_Tc.mass                # The mass of the 'representing' isotope (in u)
el_U.weight               # Abundance-weighted average mass of the stable isotopes
                          # If unstable: mass of the most stable isotope (in u)
el_C.representing_isotope # The most abundant or most stable isotope
el_At.is_stable           # Returns True if one isotope of this element is
                          # stable [Also True for Th, Pa, U]
```
**Additional properties** of an element `elem = PT.element(z)` (where `z` is an
 atomic number between 1 and 118) can be accessed via:

```python
elem.density              # in g/cm³
elem.melting_point        # at standard pressure (in K)
elem.boiling_point        # at standard pressure (in K)
elem.electro_negativity   # EN according to the Pauling scale
elem.abundance_crust      # elemental abundance in the earth's crust (in mg/kg)
elem.covalent_r_single    # Pyykkö's covalent radius for a single bond (in pm)
elem.covalent_r_double    # Pyykkö's covalent radius for a double bond (in pm)
elem.covalent_r_triple    # Pyykkö's covalent radius for a triple bond (in pm)
elem.vdW_r                # van-der-Waals radius (in pm)
elem.e_config             # html string for the electron configuration
                          # (e.g. "1s<sup>1</sup>" for H)
```

*Note that not all properties are available for every element.*
*If a property is unavailable, its value is set to 0.*

We can also access isotopes through the `PT.element().isotopes()` method
 returns a list of all isotopes or through elem.representing_isotope (returning
 only one representative).

```python
iso_D = PT.element(1).isotope(2)

# This is the same as PT.element(6).isotope(12)
iso_C12 = el_C.representing_isotope

iso_F = PT.element('F').isotope(19)
iso_Tc = PT.element(43).isotope(100)
iso_U = PT.element('U').isotope(235)
```

Internally this is handled through an `Nuclide` class which has the following
 properties:

```python
iso_D.atomic_symbol       # Returns 'D' in this case
iso_C12.atomic_number     # Returns 6 in this case
iso_F.mass_number         # Number of nucleons (aka mass number), 19 in this case
iso_Tc.atomic_mass        # mass of the nuclide
iso_U.abundance           # Nat. Isotopic abundance (0 < a <= 1)
                          # (if this is exactly zero it is unstable)
```


### Constants


In your python import section put:

```python
from chemphysconst import Constants
```

Initialise:

```python
CONST = Constants()
```
The Constants library provides access to often used physical constants.
Currently the following are implemented.

```python
CONST.planck_constant(unit)      # Plack constant (unit = "J*s" | "eV*s" | "Eh*s" )
CONST.electron_mass(unit)        # electron mass (unit = "u" | "kg")
CONST.neutron_mass(unit)         # neutron mass (unit = "u" | "kg")
CONST.proton_mass(unit)          # proton mass (unit = "u" | "kg")
CONST.avogadro_constant()        # Avogadro constant in mol^-1
CONST.bohr_magneton()            # Bohr magneton in J/T.
CONST.bohr_radius(unit)          # Bohr radius (unit = "m" | "Angstroem")
CONST.Boltzmann_constant(unit)   # Boltzmann constant (unit = "J/K" | "eV/K" | "cm^-1/K")
CONST.hartree_energy(unit)       # Hartree energy (unit = "kJ/mol" | "J" | "kcal/mol" | "eV")
CONST.atomic_mass_constant(unit) # atomic mass unit (unit = 'kg' | 'atomic')
CONST.electron_volt()            # eV in J
CONST.speed_of_light()           # speed of light in vacuum in m/s
CONST.standard_atmosphere()      # standard atmosphere in Pa
CONST.standard_temperature()     # standard temperature in K
CONST.molar_gas_constant()       # molar gas constant R in J mol^-1 K^-1
```

If there are multiple units to choose from, the first one is the default
 option.
All of the above mentioned return values are of numpy.float type.

Special mathematical definitions currently only hold the
 (3rd order Levi-Civita tensor)[https://en.wikipedia.org/wiki/Levi-Civita_symbol#Three_dimensions]
 as a numpy array:

```python
CONST.third_order_LeviCevita_Tensor()
```


## Sources

This project came about when working with the
 [PeriodicTable package](http://www.reflectometry.org/danse/elements.html).
The data is sourced from a
 [NIST compilation](http://www.nist.gov/pml/data/comp.cfm)
 compiling atomic weights and isotopic compositions for all elements
 (state 2015).
Additional properties were originally sourced from the
 [List of elements](https://en.wikipedia.org/wiki/List_of_elements)
 Wikipedia page, but are now self-curated.
Pyykkö`s covalent radii are obtained from three of his papers:
 [Single-Bond Covalent Radii](https://dx.doi.org/10.1002/chem.200800987),
 [Double-Bond Covalent Radii](https://dx.doi.org/10.1002/chem.200901472),
 [Triple-Bond Covalent Radii](https://dx.doi.org/10.1002/chem.200401299).
The van-der-radii are taken from  CRC Handbook of Chemistry and Physics,
*96th ed.*, Section:'Atomic Radii of the Elements' (except for Hydrogen, where
 120 pm instead of 110 pm is used).

All data for the physical constants is obtained from the
 [complete listing](http://physics.nist.gov/cuu/Constants/Table/allascii.txt) of
 the **CODATA internationally recommended fundamental physical constants 2014**
 curated by (NIST)[http://physics.nist.gov/cuu/Constants/index.html]
