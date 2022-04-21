import unittest

from common.evidence import get_exponent, get_mantissa


class TestExponent(unittest.TestCase):
    def test_exponent_positive(self):
        self.assertEqual(get_exponent(14.2857), 1, "Should be 1")

    def test_exponent_negative(self):
        self.assertEqual(get_exponent(-1.26e-22), -22, "Should be -22")
        self.assertEqual(get_exponent(8.06e-08), -8, "Should be -8")
        self.assertEqual(get_exponent(0.142857), -1, "Should be -1")


class TestMantissa(unittest.TestCase):
    def test_mantissa_positive(self):
        self.assertEqual(get_mantissa(14.2857), 1.429, "Should be 1.429")
        self.assertEqual(get_mantissa(0.142857), 1.429, "Should be 1.429")
        self.assertEqual(get_mantissa(8.06e-08), 8.06, "Should be 8.06")

    def test_mantiss_negative(self):
        self.assertEqual(get_mantissa(-1.26e-22), -1.26, "Should be -1.26")

    def test_mantissa_large(self):
        self.assertEqual(get_mantissa(7.364654646546544646e-44), 7.365, "Should be 7.365")


if __name__ == '__main__':
    unittest.main()
