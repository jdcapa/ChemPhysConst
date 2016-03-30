import unittest

from periodic_table.nuclides import PeriodicTable


class NuclideTest(unittest.TestCase):

    def dict_contains(self, dictionary, low=1, high=92):
        '''
        Tests if all numbers from low to high are contained in a dictionary.
        '''
        for i in range(low, high + 1):
            if i not in dictionary:
                return False
        return True

    def test_atomic_number_to_symbol(self):
        PT = PeriodicTable()
        an_to_sy = PT.atomic_number_to_symbol()
        self.assertEqual(an_to_sy[1], "H")
        self.assertEqual(an_to_sy[12], "Mg")
        self.assertEqual(an_to_sy[77], "Ir")
        self.assertEqual(an_to_sy[92], "U")
        self.assertEqual(an_to_sy[116], "Lv")
        self.assertTrue(self.dict_contains(an_to_sy))

    def test_symbol_to_atomic_number(self):
        PT = PeriodicTable()
        sy_to_an = PT.symbol_to_atomic_number()
        self.assertEqual(sy_to_an["H"], 1)
        self.assertEqual(sy_to_an["Mg"], 12)
        self.assertEqual(sy_to_an["Ir"], 77)
        self.assertEqual(sy_to_an["U"], 92)
        self.assertEqual(sy_to_an["Lv"], 116)

    def test_element_isotopic_data(self):
        PT = PeriodicTable()
        # Deuterium
        iso_H = PT.element_isotopic_data(1)
        self.assertEqual(iso_H[2]['atomic number'], 1)
        self.assertEqual(iso_H[2]['atomic symbol'], 'D')
        self.assertEqual(iso_H[2]['isotopic abundance'], 0.000115)
        self.assertEqual(iso_H[2]['atomic mass'], 2.01410177812)
        # self.assertEqual(iso_H['atomic weight'], 1.007975)
        # self.assertTrue(iso_H['stable'])
        # Fluorine-19
        iso_F = PT.element_isotopic_data(9)
        self.assertEqual(iso_F[19]['atomic number'], 9)
        self.assertEqual(iso_F[19]['atomic symbol'], 'F')
        self.assertEqual(iso_F[19]['isotopic abundance'], 1.0)
        self.assertEqual(iso_F[19]['atomic mass'], 18.99840316273)
        # self.assertEqual(iso_F['atomic weight'], 18.998403163)
        # self.assertTrue(iso_F['stable'])
        # Technetium-100
        isoTc = PT.element_isotopic_data(43)
        self.assertEqual(isoTc[100]['atomic number'], 43)
        self.assertEqual(isoTc[100]['atomic symbol'], 'Tc')
        self.assertEqual(isoTc[100]['isotopic abundance'], 0.0)
        self.assertEqual(isoTc[100]['atomic mass'], 99.9076539)
        # self.assertEqual(isoTc['atomic weight'], 97.9072124)
        # self.assertFalse(isoTc['stable'])
        # Uranium-235
        iso_U = PT.element_isotopic_data(92)
        self.assertEqual(iso_U[235]['atomic number'], 92)
        self.assertEqual(iso_U[235]['atomic symbol'], 'U')
        self.assertEqual(iso_U[235]['isotopic abundance'], 0.007204)
        self.assertEqual(iso_U[235]['atomic mass'], 235.0439301)
        # self.assertEqual(iso_U['atomic weight'], 238.02891)
        # self.assertTrue(iso_U['stable'])


if __name__ == '__main__':
    unittest.main()
