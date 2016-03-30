PeriodicTable
======================

A small python library allowing access to atomic weights and isotopic
 compositions for all known chemical elements.

Installation *(not yet working)*
------------

**Requirements:** Python 2.7, pip

Simply use pip to install the periodic table python library:

    $ pip install https://github.com/jdcapa/PeriodicTable

Usage
-----

In your python import section put:

`from periodic_table.nuclides import PeriodicTable`

Initialise:

`PT = PeriodicTable()`

Access elements via:

`elements = PT.get_elements()`

The elements dictionary contains loads of properties for the elements
 themselves:

```python
elements['H']['atomic number'] 
elements['He']['representing isotope']  # The most abundant or most stable isotope  
elements['C']['representing mass']  # The mass of the 'representing' isotope
elements['At']['is stable']  # Returns True if one isotope of this element is stable 
                             #  [Also True for Th, Pa, U]
elements['Cl']['atomic symbol'] 
elements['U']['atomic weight']  # Abundance-weighted average mass of the stable isotopes
```
as well as their isotopes:

```python
iso_Cl = elements['Cl']['isotopes']  # A dictionary of all isotopes
iso_Cl[35]['atomic mass']  # Atomic mass of a certain isotopes
iso_Cl[37]['isotopic abundance']  # Isotopic abundance 
                                  # (if this is exactly zero it is unstable)
```

Further functions include ``PT.atomic_number_to_symbol()`` and 
 ``PT.symbol_to_atomic_number()`` which return mapping dictionaries for 
 atomic numbers <=> atomic symbols.


Sources
-------

This project came about when working with the
[PeriodicTable package](http://www.reflectometry.org/danse/elements.html).
The data is sourced from a
[NIST compilation](http://www.nist.gov/pml/data/comp.cfm)
compiling atomic weights and isotopic compositions for all elements
(state 2015).

