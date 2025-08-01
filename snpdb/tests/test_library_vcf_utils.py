import tempfile

import cyvcf2
from django.conf import settings
from django.test import TestCase

from library.genomics.vcf_utils import write_vcf_from_variant_coordinates, vcf_to_variant_coordinates, \
    vcf_to_variant_coordinates_and_records, vcf_get_ref_alt_svlen_and_modification
from snpdb.models import VariantCoordinate, Sequence
from upload.models import ModifiedImportedVariant


class TestVCFUtils(TestCase):
    """ Test for library.utils - Django unit tests are run in apps so this needs to be here """
    def test_write_vcf_no_ids(self):
        variant_coordinates = [
            VariantCoordinate(chrom='1', position=123456, ref='A', alt='T'),
            VariantCoordinate(chrom='1', position=123457, ref='G', alt='T'),
            VariantCoordinate(chrom='2', position=123456, ref='A', alt='T'),
            VariantCoordinate(chrom='3', position=123456, ref='A', alt='T'),
        ]

        with tempfile.NamedTemporaryFile(delete=True) as temp_file:
            write_vcf_from_variant_coordinates(temp_file.name, variant_coordinates)
            out_vc = list(vcf_to_variant_coordinates(temp_file.name))
            self.assertEqual(variant_coordinates, out_vc)

    def test_write_vcf_with_ids(self):
        variant_coordinates = [
            VariantCoordinate(chrom='1', position=123456, ref='A', alt='T'),
            VariantCoordinate(chrom='1', position=123457, ref='G', alt='T'),
            VariantCoordinate(chrom='2', position=456789, ref='A', alt='T'),
            VariantCoordinate(chrom='3', position=789456, ref='A', alt='T'),
        ]
        vcf_ids = [
            "1_123456",
            "1_123457",
            "2_456789",
            "3_789456",
        ]
        with tempfile.NamedTemporaryFile(delete=True) as temp_file:
            write_vcf_from_variant_coordinates(temp_file.name, variant_coordinates=variant_coordinates, vcf_ids=vcf_ids)
            out_vc_r = vcf_to_variant_coordinates_and_records(temp_file.name)
            for (in_vc, in_id), (out_vc, out_record) in zip(zip(variant_coordinates, vcf_ids), out_vc_r):
                self.assertEqual(in_vc, out_vc)
                self.assertEqual(in_id, out_record.ID)

    def test_write_vcf_with_integer_ids(self):
        variant_coordinates = [
            VariantCoordinate(chrom='1', position=123456, ref='A', alt='T'),
        ]
        vcf_ids = [
            42,
        ]
        with tempfile.NamedTemporaryFile(delete=True) as temp_file:
            write_vcf_from_variant_coordinates(temp_file.name, variant_coordinates=variant_coordinates,
                                               vcf_ids=vcf_ids)
            vcf_to_variant_coordinates_and_records(temp_file.name)

    def test_write_vcf_sort_coordinates(self):
        """ Make sure that if you have out of order coordinates they are sorted, and IDs still match """
        variant_coordinates = [
            VariantCoordinate(chrom='1', position=123457, ref='G', alt='T'),
            VariantCoordinate(chrom='1', position=123456, ref='A', alt='T'),
            VariantCoordinate(chrom='3', position=789456, ref='A', alt='T'),
            VariantCoordinate(chrom='2', position=456789, ref='A', alt='T'),
        ]
        vcf_ids = [
            "1_123457",
            "1_123456",
            "3_789456",
            "2_456789",
        ]
        with tempfile.NamedTemporaryFile(delete=True) as temp_file:
            write_vcf_from_variant_coordinates(temp_file.name, variant_coordinates=variant_coordinates, vcf_ids=vcf_ids)
            previous_out_vc = None
            for out_vc, out_record in vcf_to_variant_coordinates_and_records(temp_file.name):
                # Check that IDs match the coordinate
                expected_id = f"{out_vc.chrom}_{out_vc.position}"
                self.assertEqual(expected_id, out_record.ID)

                # Check that they are in order
                if previous_out_vc is not None:
                    self.assertLess(previous_out_vc, out_vc)
                previous_out_vc = out_vc

    def test_write_vcf_has_end(self):
        variant_coordinates = [
            VariantCoordinate(chrom='1', position=123456, ref='A', alt='<DEL>', svlen=5000),
        ]

        with tempfile.NamedTemporaryFile(delete=True) as temp_file:
            write_vcf_from_variant_coordinates(temp_file.name, variant_coordinates)
            out = vcf_to_variant_coordinates_and_records(temp_file.name)
            for in_vc, (out_vc, out_record) in zip(variant_coordinates, out):
                self.assertEqual(out_record.INFO["END"], in_vc.end)

    def test_vcf_split_multi_get_ref_alt_svlen_and_modification(self):
        filename = "snpdb/tests/test_data/svlen_split_multi_allele.vcf"
        for v in cyvcf2.Reader(filename):
            _ref, _alt, svlen, _modification = vcf_get_ref_alt_svlen_and_modification(v, old_variant_info=ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG)
            self.assertTrue(isinstance(svlen, int))

    def _assert_del_svlen(self, alt, svlen):
        if Sequence.allele_is_symbolic(alt):
            self.assertTrue(isinstance(svlen, int))
            if alt == "<DEL>":
                if settings.VARIANT_SYMBOLIC_ALT_SVLEN_ALWAYS_POSITIVE:
                    self.assertGreater(svlen, 0)
                else:
                    self.assertLess(svlen, 0)

    def test_vcf_manta_get_ref_alt_svlen_and_modification(self):
        filename = "snpdb/tests/test_data/manta.vcf"
        for v in cyvcf2.Reader(filename):
            _ref, alt, svlen, _modification = vcf_get_ref_alt_svlen_and_modification(v, old_variant_info=ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG)
            self._assert_del_svlen(alt, svlen)

    def test_vcf_no_svlen_del_only_end_get_ref_alt_svlen_and_modification(self):
        filename = "snpdb/tests/test_data/symbolic_alt_del_no_svlen_only_end.vcf"
        for v in cyvcf2.Reader(filename):
            _ref, alt, svlen, _modification = vcf_get_ref_alt_svlen_and_modification(v, old_variant_info=ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG)
            self._assert_del_svlen(alt, svlen)

