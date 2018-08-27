"""Contains tests and tests suits for all methods in DivEn project
"""

import sys
import os
import unittest
import io
import time

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

    def test_zero_approximation(self):
        """Testing zero order approximation
        """

        ans = [(0, 0, 0, 0), (0, 0, 1, 3), (0, 1, 0, 2), (0, 2, 0, 4), (1, 0, 0, 1), (1, 0, 1, 4),
               (1, 1, 0, 3), (2, 0, 0, 2), (2, 1, 0, 4), (3, 0, 0, 3), (4, 0, 0, 4)]
        self.assertListEqual(list(RSPT.gen_zero_approximation(0, 4, [1, 2, 3])), ans)
        states = RSPT.zero_approximation(200, [1, 2, 3])
        self.assertEqual(len(states), 223872)


class TestHarmOscillatorMethods(unittest.TestCase):
    """Class for testing Harmonic oscillator methods
    """

    def test_zero_elements(self):
        """Test for correctness for defining zero elements
        """

        self.assertEqual(harm_oscill.ZeroEl(0, 0), False)
        self.assertEqual(harm_oscill.ZeroEl(2, 0), False)
        self.assertEqual(harm_oscill.ZeroEl(4, 2), False)
        self.assertEqual(harm_oscill.ZeroEl(1, 1), False)
        self.assertEqual(harm_oscill.ZeroEl(3, 1), False)
        self.assertEqual(harm_oscill.ZeroEl(3, 3), False)
        self.assertEqual(harm_oscill.ZeroEl(8, 4), False)
        self.assertEqual(harm_oscill.ZeroEl(2, 3), True)
        self.assertEqual(harm_oscill.ZeroEl(2, 4), True)
        self.assertEqual(harm_oscill.ZeroEl(4, 3), True)
        self.assertEqual(harm_oscill.ZeroEl(2, 1), True)
        self.assertEqual(harm_oscill.ZeroEl(9, 3), True)


if __name__ == '__main__':
    unittest.main(verbosity=2)
