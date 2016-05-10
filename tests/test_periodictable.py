import unittest

from chemphysconst import PeriodicTable


class PeriodicTableTest(unittest.TestCase):

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
        iso_H = PT.element(1).isotope(2)
        iso_F = PT.element(9).isotope(19)
        iso_Tc = PT.element(43).isotope(100)
        iso_U = PT.element(92).isotope(235)

        # Deuterium
        self.assertEqual(iso_H.atomic_number, 1)
        self.assertEqual(iso_H.atomic_symbol, 'D')
        self.assertEqual(iso_H.abundance, 0.000115)
        self.assertEqual(iso_H.atomic_mass, 2.01410177812)
        # Fluorine-19
        self.assertEqual(iso_F.atomic_number, 9)
        self.assertEqual(iso_F.atomic_symbol, 'F')
        self.assertEqual(iso_F.abundance, 1.0)
        self.assertEqual(iso_F.atomic_mass, 18.99840316273)
        # Technetium-100
        self.assertEqual(iso_Tc.atomic_number, 43)
        self.assertEqual(iso_Tc.atomic_symbol, 'Tc')
        self.assertEqual(iso_Tc.abundance, 0.0)
        self.assertEqual(iso_Tc.atomic_mass, 99.9076539)
        # Uranium-235
        self.assertEqual(iso_U.atomic_number, 92)
        self.assertEqual(iso_U.atomic_symbol, 'U')
        self.assertEqual(iso_U.abundance, 0.007204)
        self.assertEqual(iso_U.atomic_mass, 235.0439301)

    def test_get_representative_isotope(self):
        PT = PeriodicTable()
        elem_H = PT.element(1)
        elem_F = PT.element(9)
        elem_Tc = PT.element(43)
        elem_U = PT.element(92)
        self.assertEqual(elem_H.representing_isotope.mass_number, 1)
        self.assertEqual(elem_F.representing_isotope.mass_number, 19)
        self.assertEqual(elem_Tc.representing_isotope.mass_number, 98)
        self.assertEqual(elem_U.representing_isotope.mass_number, 238)

    def test_stability(self):
        PT = PeriodicTable()
        elem_H = PT.element(1)
        elem_F = PT.element(9)
        elem_Tc = PT.element(43)
        elem_U = PT.element(92)
        self.assertTrue(elem_H.is_stable)
        self.assertTrue(elem_F.is_stable)
        self.assertFalse(elem_Tc.is_stable)
        self.assertTrue(elem_U.is_stable)

    def test_calculate_atomic_weight(self):
        PT = PeriodicTable()
        elem_H = PT.element(1)
        elem_F = PT.element(9)
        elem_Tc = PT.element(43)
        elem_U = PT.element(92)

        self.assertEqual(elem_H.weight, 1.0079407540557772)
        self.assertEqual(elem_F.weight, 18.99840316273)
        self.assertEqual(elem_Tc.weight, 97.9072124)
        self.assertEqual(elem_U.weight, 238.0289104616574)


if __name__ == '__main__':
    unittest.main()
