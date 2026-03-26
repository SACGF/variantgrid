"""
Tests for pedigree.ped.export_ped and round-trips through import_ped.

Two complementary strategies:
- File-content tests (TestWriteTrioPedFileContents): check the exact PED encoding written to disk.
  These catch bugs that a round-trip cannot — e.g. if both export and import have the same
  wrong encoding, the round-trip passes but Somalier (which reads the raw file) still gets
  wrong phenotype data.
- Round-trip tests (TestWriteTrioPedRoundTrip): check end-to-end behaviour after import.
  Used for properties where the encoding is ambiguous ('M'/'F' vs '1'/'2' both import fine).
"""
import os
import tempfile

import django.test
from django.contrib.auth.models import User

from patients.models_enums import Sex
from pedigree.ped.export_ped import write_trio_ped, write_unrelated_ped
from pedigree.ped.import_ped import import_ped


def _import_file(path, user):
    with open(path) as f:
        return import_ped(f, "test.ped", user)


class _TmpPedFile(django.test.TestCase):
    """Mixin: creates a temp file for each test and deletes it in tearDown."""

    def setUp(self):
        tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.ped', delete=False)
        tmp.close()
        self.fname = tmp.name

    def tearDown(self):
        os.unlink(self.fname)

    def _read_rows(self):
        """Return dict of sample_name → list of PED columns."""
        with open(self.fname) as f:
            rows = [line.strip().split('\t') for line in f if line.strip()]
        return {row[1]: row for row in rows}


class TestWriteTrioPedFileContents(_TmpPedFile):
    """Check the raw affection encoding in the written PED file.
    Standard PED: 2=affected, 1=unaffected, 0=unknown.
    BUG-1: write_trio_ped currently uses 1=affected, 0=unaffected — opposite of the standard.
    """

    def _write_default(self, father_affected=False, mother_affected=False):
        write_trio_ped(self.fname, 'proband', Sex.MALE,
                       'father', father_affected, 'mother', mother_affected,
                       family_code='FAM001')

    def test_proband_affection_column_is_two(self):
        # Proband is always affected; must be '2' in the PED file
        self._write_default()
        rows = self._read_rows()
        self.assertEqual(rows['proband'][5], '2',
                         "Proband affection should be '2' (affected) — currently writes '1' (bug)")

    def test_affected_parent_affection_column_is_two(self):
        # father_affected=True must write '2', not '1'
        self._write_default(father_affected=True)
        rows = self._read_rows()
        self.assertEqual(rows['father'][5], '2',
                         "Affected father affection should be '2' — currently writes '1' (bug)")

    def test_unaffected_parent_affection_column_is_one(self):
        # father_affected=False must write '1' (unaffected), not '0' (unknown)
        self._write_default()
        rows = self._read_rows()
        self.assertEqual(rows['father'][5], '1',
                         "Unaffected father affection should be '1' — currently writes '0' (bug)")


class TestWriteTrioPedRoundTrip(_TmpPedFile):
    """End-to-end: write_trio_ped → import_ped → verify model state.
    Used for sex and parent links where the exact encoding format doesn't matter
    as long as it round-trips correctly.
    """

    def setUp(self):
        super().setUp()
        self.user = User.objects.create_user('export_test')

    def _write_and_import(self, proband_sex=Sex.MALE, father_affected=False, mother_affected=False):
        write_trio_ped(self.fname, 'proband', proband_sex,
                       'father', father_affected, 'mother', mother_affected,
                       family_code='FAM001')
        return _import_file(self.fname, self.user)

    def test_round_trip_passes_validation(self):
        # BUG-1 manifests here: wrong affection → no affected individual → ValidationError
        _, families = self._write_and_import()
        self.assertEqual(families[0].errors, [])

    def test_proband_is_affected_after_round_trip(self):
        _, families = self._write_and_import()
        proband = families[0].pedfilerecord_set.get(sample='proband')
        self.assertTrue(proband.affection)

    def test_sex_round_trip(self):
        _, families = self._write_and_import(proband_sex=Sex.MALE)
        records = {r.sample: r for r in families[0].pedfilerecord_set.all()}
        self.assertEqual(records['proband'].sex, Sex.MALE)
        self.assertEqual(records['father'].sex, Sex.MALE)
        self.assertEqual(records['mother'].sex, Sex.FEMALE)

    def test_parent_links_round_trip(self):
        _, families = self._write_and_import()
        proband = families[0].pedfilerecord_set.get(sample='proband')
        self.assertIsNotNone(proband.father)
        self.assertEqual(proband.father.sample, 'father')
        self.assertIsNotNone(proband.mother)
        self.assertEqual(proband.mother.sample, 'mother')


class TestWriteUnrelatedPed(_TmpPedFile):

    def test_creates_one_row_per_sample(self):
        write_unrelated_ped(self.fname, ['s1', 's2', 's3'], family_code='FAM001')
        rows = self._read_rows()
        self.assertEqual(set(rows.keys()), {'s1', 's2', 's3'})

    def test_no_parent_columns(self):
        # Unrelated individuals must have '0' for both father and mother columns
        write_unrelated_ped(self.fname, ['s1'], family_code='FAM001')
        rows = self._read_rows()
        self.assertEqual(rows['s1'][2], '0')
        self.assertEqual(rows['s1'][3], '0')

    def test_empty_sample_list_creates_empty_file(self):
        write_unrelated_ped(self.fname, [], family_code='FAM001')
        rows = self._read_rows()
        self.assertEqual(rows, {})
