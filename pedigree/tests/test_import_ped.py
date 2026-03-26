"""
Tests for pedigree.ped.import_ped — import pipeline (DB required).
"""
import io

import django.test
from django.contrib.auth.models import User

from patients.models_enums import Sex
from pedigree.models import PedFile
from pedigree.ped.import_ped import import_ped
from snpdb.models import ImportStatus


def _import(ped_text, user, name="test.ped"):
    return import_ped(io.StringIO(ped_text), name, user)


VALID_TRIO = (
    "FAM001 proband father mother 1 2\n"
    "FAM001 father  0      0      1 1\n"
    "FAM001 mother  0      0      2 1\n"
)


class TestImportPedBasic(django.test.TestCase):

    def setUp(self):
        self.user = User.objects.create_user('import_test')

    def test_valid_trio_imports_successfully(self):
        ped_file, families = _import(VALID_TRIO, self.user)
        self.assertEqual(ped_file.import_status, ImportStatus.SUCCESS)
        self.assertEqual(len(families), 1)

    def test_proband_is_affected_after_import(self):
        _, families = _import(VALID_TRIO, self.user)
        proband = families[0].pedfilerecord_set.get(sample='proband')
        self.assertTrue(proband.affection)

    def test_parent_relationships_wired_correctly(self):
        _, families = _import(VALID_TRIO, self.user)
        proband = families[0].pedfilerecord_set.get(sample='proband')
        self.assertIsNotNone(proband.father)
        self.assertEqual(proband.father.sample, 'father')
        self.assertIsNotNone(proband.mother)
        self.assertEqual(proband.mother.sample, 'mother')

    def test_sex_parsed_correctly(self):
        _, families = _import(VALID_TRIO, self.user)
        records = {r.sample: r for r in families[0].pedfilerecord_set.all()}
        self.assertEqual(records['proband'].sex, Sex.MALE)
        self.assertEqual(records['father'].sex, Sex.MALE)
        self.assertEqual(records['mother'].sex, Sex.FEMALE)

    def test_two_family_ped_creates_two_families(self):
        ped_text = (
            "FAM001 proband_a 0 0 1 2\n"
            "FAM002 proband_b 0 0 2 2\n"
        )
        ped_file, families = _import(ped_text, self.user)
        self.assertEqual(len(families), 2)
        self.assertEqual(ped_file.import_status, ImportStatus.SUCCESS)


class TestImportPedValidation(django.test.TestCase):

    def setUp(self):
        self.user = User.objects.create_user('import_validation_test')

    def test_no_affected_raises(self):
        ped_text = (
            "FAM001 proband 0 0 1 1\n"
            "FAM001 father  0 0 1 1\n"
        )
        with self.assertRaises(Exception):
            _import(ped_text, self.user)

    def test_no_affected_sets_error_status(self):
        # Distinct from test_no_affected_raises: checks the DB is updated, not just the exception
        ped_text = (
            "FAM001 proband 0 0 1 1\n"
            "FAM001 father  0 0 1 1\n"
        )
        try:
            _import(ped_text, self.user)
        except Exception:
            pass
        ped_file = PedFile.objects.filter(name="test.ped", user=self.user).last()
        self.assertEqual(ped_file.import_status, ImportStatus.ERROR)

    def test_wrong_sex_father_raises(self):
        ped_text = (
            "FAM001 proband father mother 1 2\n"
            "FAM001 father  0      0      2 1\n"  # sex=2 (Female)
            "FAM001 mother  0      0      2 1\n"
        )
        with self.assertRaises(Exception):
            _import(ped_text, self.user)

    def test_phenotips_minus9_affection_imports_as_none(self):
        # '-9' is a Phenotips-specific unknown affection code — must not raise
        ped_text = (
            "FAM001 proband 0 0 1 2\n"
            "FAM001 sibling 0 0 1 -9\n"
        )
        _, families = _import(ped_text, self.user)
        sibling = families[0].pedfilerecord_set.get(sample='sibling')
        self.assertIsNone(sibling.affection)

    def test_missing_parent_reference_raises(self):
        # Parent ID present in a record but not defined as its own row
        ped_text = "FAM001 proband ghost_dad 0 1 2\n"
        with self.assertRaises(Exception):
            _import(ped_text, self.user)
