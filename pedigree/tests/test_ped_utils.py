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

    def test_unknown_value_returns_none_not_raises(self):
        # Asymmetric contract: get_sex silently returns None while get_affection raises
        # A PED file with an unrecognised sex value is silently treated as unknown sex
        self.assertIsNone(get_sex('99'))

    def test_lowercase_m_not_recognised(self):
        # BUG-4: 'm' is not in SEX_LOOKUP — silently returns None instead of MALE
        result = get_sex('m')
        self.assertEqual(result, Sex.MALE,
                         "Lowercase 'm' should map to MALE — currently returns None (bug)")

    def test_lowercase_f_not_recognised(self):
        # BUG-4: 'f' is not in SEX_LOOKUP — silently returns None instead of FEMALE
        result = get_sex('f')
        self.assertEqual(result, Sex.FEMALE,
                         "Lowercase 'f' should map to FEMALE — currently returns None (bug)")


class TestGetParentId(unittest.TestCase):

    def test_integer_zero_is_unknown(self):
        # pandas may supply int 0 (not string '0') for an all-numeric parent column
        self.assertIsNone(get_parent_id(0))

    def test_dot_is_unknown(self):
        # '.' is a less obvious alternative unknown-parent marker used by some tools
        self.assertIsNone(get_parent_id('.'))

    def test_nan_is_unknown(self):
        # BUG-6: float('nan') not in UNKNOWN_PARENT_VALUES — passes through as a "valid" ID
        result = get_parent_id(float('nan'))
        self.assertIsNone(result,
                          "NaN should be treated as an unknown parent, not a valid ID")
