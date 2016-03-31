# Just so one can update the HDF5 database if the NIST data is changed.
from periodic_table.periodictable import ExportPeriodicTableNIST

EPTN = ExportPeriodicTableNIST()
EPTN.export_elements_to_hdf5()