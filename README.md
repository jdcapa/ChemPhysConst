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

`elements['H']['atomic number']`

`elements['He']['representing isotope']`

`elements['C']['representing mass']`

`elements['At']['is stable']  # Returns True if`

`elements['Cl']['atomic symbol']`

`elements['U']['atomic weight']`




Sources
-------

This project came about when working with the
[PeriodicTable package](http://www.reflectometry.org/danse/elements.html).
The data is sourced from a
[NIST compilation](http://www.nist.gov/pml/data/comp.cfm)
compiling atomic weights and isotopic compositions for all elements
(state 2015).

