PeriodicTable
======================

A small python library allowing access to atomic weights and isotopic
 compositions for all known chemical elements.
It is meant to be lightweight and fast due to the hdf5 database. 

Installation *(not yet working)*
------------

**Requirements:** Python 2.7, pip,
 [hd5py](http://docs.h5py.org/en/latest/index.html)

Simply use pip to install the periodic table python library:

    $ pip install https://github.com/jdcapa/PeriodicTable

Usage
-----

In your python import section put:

```python
from periodic_table.periodictable import PeriodicTable
```

Initialise:

```python
PT = PeriodicTable()
```

where PT contains the dictionaries

```python
PT.an_to_sy  # mapping all the atomic numbers to the symbols
PT.sy_to_an  # mapping all the symbols to the atomic numbers
```

Access elements via the PT.elements() method:

```python
el_H = PT.element(1)    # Either via the atomic number (Hydrogen)
el_F = PT.element('F')  # or the atomic symbol (Fluorine)
el_Tc = PT.element(43)  # Technetium
el_U = PT.element('U')  # Uranium
```

The elements contains loads of properties themselves:

```python
el_H.symbol                 # atomic number
el_F.number                 # atomic symbol
el_Tc.mass                  # The mass of the 'representing' isotope
el_U.weight                 # Abundance-weighted average mass of the stable isotopes
                            #  If unstable: mass of the most stable isotope
el_C.representing_isotope   # The most abundant or most stable isotope
el_At.is_stable             # Returns True if one isotope of this element is
                            #  stable [Also True for Th, Pa, U]
```
as well as their isotopes (through the PT.elements().isotopes() method):

```python
iso_D = PT.element(1).isotope(2)
iso_F = PT.element('F').isotope(19)
iso_Tc = PT.element(43).isotope(100)
iso_U = PT.element('U').isotope(235)
```

which in turn have properties:

```python
iso_D.atomic_symbol         # Returns 'D' in this case
iso_F.atomic_number
iso_F.mass_number           # Number of nucleons (aka *mass number*)
iso_Tc.atomic_mass          # mass of the nuclide
iso_U.abundance             # Nat. Isotopic abundance (0 < a <= 1)
                            # (if this is exactly zero it is unstable)
```


Sources
-------

This project came about when working with the
[PeriodicTable package](http://www.reflectometry.org/danse/elements.html).
The data is sourced from a
[NIST compilation](http://www.nist.gov/pml/data/comp.cfm)
compiling atomic weights and isotopic compositions for all elements
(state 2015).

