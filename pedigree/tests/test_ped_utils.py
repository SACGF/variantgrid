"""
Tests for pedigree.ped.ped_file_utils — pure functions, no DB required.
"""
import unittest

from patients.models_enums import Sex
from pedigree.ped.ped_file_utils import get_affection, get_parent_id, get_sex


class TestGetAffection(unittest.TestCase):

    def test_minus_nine_phenotips_is_unknown(self):
        # Non-obvious: Phenotips uses '-9' for unknown, not '0'
        self.assertIsNone(get_affection('-9'))

    def test_unknown_value_raises(self):
        # Strict contract: unrecognised values raise rather than silently returning None
        with self.assertRaises(ValueError):
            get_affection('99')


class TestGetSex(unittest.TestCase):
    # Asymmetric contract: get_sex silently returns None while get_affection raises
    # A PED file with an unrecognised sex value is silently treated as unknown sex

    def test_lowercase_m_is_male(self):
        self.assertEqual(get_sex('m'), Sex.MALE)

    def test_lowercase_f_is_female(self):
        self.assertEqual(get_sex('f'), Sex.FEMALE)


class TestGetParentId(unittest.TestCase):

    def test_integer_zero_is_unknown(self):
        # pandas may supply int 0 (not string '0') for an all-numeric parent column
        self.assertIsNone(get_parent_id(0))

    def test_dot_is_unknown(self):
        # '.' is a less obvious alternative unknown-parent marker used by some tools
        self.assertIsNone(get_parent_id('.'))

    def test_nan_is_unknown(self):
        # pandas supplies NaN for a blank cell in an all-numeric parent column
        self.assertIsNone(get_parent_id(float('nan')),
                          "NaN should be treated as an unknown parent, not a valid ID")
