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
        iso_H = PT.element_isotopic_data(1)
        iso_F = PT.element_isotopic_data(9)
        iso_Tc = PT.element_isotopic_data(43)
        iso_U = PT.element_isotopic_data(92)

        # Deuterium
        self.assertEqual(iso_H[2]['atomic number'], 1)
        self.assertEqual(iso_H[2]['atomic symbol'], 'D')
        self.assertEqual(iso_H[2]['isotopic abundance'], 0.000115)
        self.assertEqual(iso_H[2]['atomic mass'], 2.01410177812)
        # Fluorine-19
        self.assertEqual(iso_F[19]['atomic number'], 9)
        self.assertEqual(iso_F[19]['atomic symbol'], 'F')
        self.assertEqual(iso_F[19]['isotopic abundance'], 1.0)
        self.assertEqual(iso_F[19]['atomic mass'], 18.99840316273)
        # Technetium-100
        self.assertEqual(iso_Tc[100]['atomic number'], 43)
        self.assertEqual(iso_Tc[100]['atomic symbol'], 'Tc')
        self.assertEqual(iso_Tc[100]['isotopic abundance'], 0.0)
        self.assertEqual(iso_Tc[100]['atomic mass'], 99.9076539)
        # Uranium-235
        self.assertEqual(iso_U[235]['atomic number'], 92)
        self.assertEqual(iso_U[235]['atomic symbol'], 'U')
        self.assertEqual(iso_U[235]['isotopic abundance'], 0.007204)
        self.assertEqual(iso_U[235]['atomic mass'], 235.0439301)

    def test_get_repr_isotope(self):
        PT = PeriodicTable()
        iso_H = PT.element_isotopic_data(1)
        iso_F = PT.element_isotopic_data(9)
        iso_Tc = PT.element_isotopic_data(43)
        iso_U = PT.element_isotopic_data(92)
        rpr_iso = PT.get_representative_isotope
        self.assertEqual(rpr_iso(iso_H), 1)
        self.assertEqual(rpr_iso(iso_F), 19)
        self.assertEqual(rpr_iso(iso_Tc), 98)
        self.assertEqual(rpr_iso(iso_U), 238)

    def test_elemtent_stability(self):
        PT = PeriodicTable()
        iso_H = PT.element_isotopic_data(1)
        iso_F = PT.element_isotopic_data(9)
        iso_Tc = PT.element_isotopic_data(43)
        iso_U = PT.element_isotopic_data(92)
        self.assertTrue(PT.elemtent_stability(iso_H))
        self.assertTrue(PT.elemtent_stability(iso_F))
        self.assertFalse(PT.elemtent_stability(iso_Tc))
        self.assertTrue(PT.elemtent_stability(iso_U))

    def test_atomic_weight(self):
        PT = PeriodicTable()
        iso_H = PT.element_isotopic_data(1)
        iso_F = PT.element_isotopic_data(9)
        iso_Tc = PT.element_isotopic_data(43)
        iso_U = PT.element_isotopic_data(92)
        aw = PT.atomic_weight
        self.assertEqual(aw(iso_H, 'range'), 1.007975)
        self.assertEqual(aw(iso_F, 'range'), 18.998403163)
        self.assertEqual(aw(iso_Tc, 'range'), 97.9072124)
        self.assertEqual(aw(iso_U, 'range'), 238.02891)

        self.assertEqual(aw(iso_H, 'average'), 1.0079407540557772)
        self.assertEqual(aw(iso_F, 'average'), 18.99840316273)
        self.assertEqual(aw(iso_Tc, 'average'), 97.9072124)
        self.assertEqual(aw(iso_U, 'average'), 238.0289104616574)

    def test_elements(self):
        PT = PeriodicTable()
        elements = PT.get_elements()
        self.assertEqual(elements['H']['atomic number'], 1)
        self.assertEqual(elements['He']['representing isotope'], 4)
        self.assertEqual(elements['C']['representing mass'], 12.0)
        self.assertFalse(elements['At']['is stable'])
        self.assertEqual(elements['Cl']['atomic symbol'], 'Cl')
        self.assertEqual(elements['U']['atomic weight'], 238.0289104616574)


if __name__ == '__main__':
    unittest.main()
