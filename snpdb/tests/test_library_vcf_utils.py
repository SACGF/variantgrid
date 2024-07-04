import tempfile

from django.test import TestCase

from library.genomics.vcf_utils import write_vcf_from_variant_coordinates, vcf_to_variant_coordinates, \
    vcf_to_variant_coordinates_and_records
from snpdb.models import VariantCoordinate


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
