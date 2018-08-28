"""Contains tests and tests suits for all methods in DivEn project
"""

import sys
import os
import unittest
import io
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../lib/RSPT'))

import supp
import mol_io
import RSPT
import harm_oscill


class TestExceptionMethods(unittest.TestCase):
    """Class for testing supplementary methods
    """

    def test_div_zero(self):
        """Testing division by zero exception
        """
        try:
            print(1/0)
        except ZeroDivisionError:
            f, n, ln, obj, err_type = supp.PrintException()
        finally:
            self.assertEqual(f, 'tests.py')
            self.assertEqual(ln, 'print(1/0)')
            self.assertEqual(str(obj), 'division by zero')
            self.assertIs(err_type, ZeroDivisionError)
            self.assertIsInstance(n, int)


    def test_index_error(self):
        """Testing index error exception (index out of range)
        """
        try:
            l = []
            print(l[0])
        except IndexError:
            f, n, ln, obj, err_type = supp.PrintException()
        finally:
            self.assertEqual(f, 'tests.py')
            self.assertEqual(ln, 'print(l[0])')
            self.assertEqual(str(obj), 'list index out of range')
            self.assertIs(err_type, IndexError)
            self.assertIsInstance(n, int)


    def test_sys_info(self):
        """Testing printing system info
        """

        f, n, ln = supp.PrintSysExitInfo()
        self.assertEqual(f, 'tests.py')
        self.assertEqual(ln, '"""')
        self.assertIsInstance(n, int)


    def test_timing_decorator(self):
        """Test timing decorator
        """

        @supp.timing
        def hold(sec):
            """Hold for 'sec' seconds
            """
            time.sleep(sec)
            return 1

        capturedOutput = io.StringIO()
        sys.stdout = capturedOutput
        hold(0.1)
        sys.stdout = sys.__stdout__
        self.assertAlmostEqual(0.1, float(capturedOutput.getvalue().split()[3]), places=1)


    def test_reader_class(self):
        """Test reader class
        """

        def g(k):
            """Test generator
            """
            for i in range(k):
                yield i
        A = supp.Reader(g(5))

        self.assertEqual(int(A.read()), 0)
        self.assertEqual(int(A.read()), 1)

        with self.assertRaises(TypeError):
            A = supp.Reader(3)

        with self.assertRaises(TypeError):
            A = supp.Reader(range(5))

        self.assertEqual(int(A.read()), 2)


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


class TestRSPTMethods(unittest.TestCase):
    """Class for testing RSPT methods
    """

    def setUp(self):
        """Setting up test objects
        """
        path = os.path.join(os.path.dirname(__file__), '../etc/NO2_test_input/')
        self.fname_freqs = path + 'Frequencies.txt'
        self.fname_coefs = path + 'Anh_coefs.txt'


    def test_zero_approximation(self):
        """Testing zero order approximation
        """

        ans = [(0, 0, 0, 0), (0, 0, 1, 3), (0, 1, 0, 2), (0, 2, 0, 4), (1, 0, 0, 1), (1, 0, 1, 4),
               (1, 1, 0, 3), (2, 0, 0, 2), (2, 1, 0, 4), (3, 0, 0, 3), (4, 0, 0, 4)]
        self.assertListEqual(list(RSPT.gen_zero_approximation(0, 4, [1, 2, 3])), ans)
        states = RSPT.zero_approximation(200, [1, 2, 3])
        self.assertEqual(len(states), 223872)


    def test_fill_perturb_mat(self):
        """Testing filling perturbation matrix
        """

        freqs = mol_io.read_freqs(self.fname_freqs)
        anh_coefs = mol_io.read_anh_coefs(self.fname_coefs)
        zero_states = RSPT.zero_approximation(10000, freqs['omega'])
        zero_states = zero_states[zero_states['v3'] % 2 == 1]
        zero_states.index = range(len(zero_states))
        print(zero_states)
        Wmat = RSPT.fill_wmat2(anh_coefs, zero_states)
        Hmat = Wmat
        Hmat.flat[::Hmat.shape[0] + 1] += zero_states['E']
        time1 = time.time()
        Eigh_values, Eigh_vectors = np.linalg.eigh(Hmat)
        time2 = time.time()
        print(Eigh_values)
        print('Matrix diagonalisation took {} s'.format(time2 - time1))



class TestHarmOscillatorMethods(unittest.TestCase):
    """Class for testing Harmonic oscillator methods
    """

    def test_mat_el_weights(self):
        """Test calculation of weights for matrix elements
        """

        Weights_dict = harm_oscill.MatElWeightDict(15)

        self.assertAlmostEqual(Weights_dict[(0, 0, 15)], 1, places=8)
        self.assertAlmostEqual(Weights_dict[(1, 1, 2)], 1, places=8)
        self.assertAlmostEqual(Weights_dict[(2, 2, 2)], 0.7071067810, places=8)
        self.assertAlmostEqual(Weights_dict[(2, 0, 2)], 2.5, places=8)
        self.assertAlmostEqual(Weights_dict[(3, 3, 2)], 0., places=8)
        self.assertAlmostEqual(Weights_dict[(3, 1, 2)], 3., places=8)
        self.assertAlmostEqual(Weights_dict[(4, 4, 2)], 0., places=8)
        self.assertAlmostEqual(Weights_dict[(4, 2, 2)], 2.121320343, places=8)
        self.assertAlmostEqual(Weights_dict[(4, 0, 2)], 9.75, places=8)
        self.assertAlmostEqual(Weights_dict[(5, 5, 9)], 21.73706511)
        self.assertAlmostEqual(Weights_dict[(5, 3, 9)], 158.7450786, places=4)
        self.assertAlmostEqual(Weights_dict[(5, 1, 9)], 432.2190198, places=4)
        self.assertAlmostEqual(Weights_dict[(6, 6, 9)], 30.74085230, places=4)
        self.assertAlmostEqual(Weights_dict[(6, 4, 9)], 309.3238593, places=4)
        self.assertAlmostEqual(Weights_dict[(6, 2, 9)], 1161.422888, places=4)
        self.assertAlmostEqual(Weights_dict[(6, 0, 9)], 2173.125000, places=4)
        self.assertAlmostEqual(Weights_dict[(7, 7, 9)], 37.64970119, places=4)
        self.assertAlmostEqual(Weights_dict[(7, 5, 9)], 532.5580950, places=4)
        self.assertAlmostEqual(Weights_dict[(7, 3, 9)], 2708.587904, places=4)
        self.assertAlmostEqual(Weights_dict[(7, 1, 9)], 6932.740046, places=4)
        self.assertAlmostEqual(Weights_dict[(8, 8, 9)], 37.64970119, places=4)
        self.assertAlmostEqual(Weights_dict[(8, 6, 9)], 799.2621590, places=4)
        self.assertAlmostEqual(Weights_dict[(8, 4, 9)], 5533.460149, places=4)
        self.assertAlmostEqual(Weights_dict[(8, 2, 9)], 18932.78406, places=4)
        self.assertAlmostEqual(Weights_dict[(8, 0, 9)], 37019.06247, places=4)


if __name__ == '__main__':
    unittest.main(verbosity=2)
