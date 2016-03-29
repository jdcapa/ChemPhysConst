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

    def test_element_map(self):
        PT = PeriodicTable()
        element_map = PT.element_map()
        self.assertEqual(element_map[1], "H")
        self.assertEqual(element_map[12], "Mg")
        self.assertEqual(element_map[77], "Ir")
        self.assertEqual(element_map[92], "U")
        self.assertEqual(element_map[116], "Lv")
        self.assertTrue(self.dict_contains(element_map))


if __name__ == '__main__':
    unittest.main()
