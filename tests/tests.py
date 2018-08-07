"""Contains tests and tests suits for all methods in DivEn project
"""

import sys
import os
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../lib/RSPT'))

import mol_io

class TestIOMethods(unittest.TestCase):
    """Class for testing IO methods
    """

    def setUp(self):
        """Setting up test objects
        """
        path = os.path.join(os.path.dirname(__file__), '../etc/HDO_test_input/')
        self.fname_freqs = path + 'Frequencies.txt'
        self.fname_coefs = path + 'Anh_coefs.txt'


    def test_basic_read(self):
        """Testing correctness of read data
        """

        freqs = mol_io.read_freqs(self.fname_freqs)
        coefs = mol_io.read_anh_coefs(self.fname_coefs)

        self.assertListEqual(freqs['omega'].tolist(), [2824.3, 1440.2, 3889.8])
        self.assertListEqual(coefs.iloc[0].tolist(), [3, 0, 0, -258.4])
        self.assertListEqual(coefs.iloc[4].tolist(), [1, 1, 1, -15.4])
        self.assertListEqual(coefs.iloc[24].tolist(), [0, 0, 4, 62])


if __name__ == '__main__':
    unittest.main(verbosity=2)
