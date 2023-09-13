from django.test import TestCase

from snpdb.models import Allele, AlleleConversionTool, GenomeBuild, Variant
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class VariantTestCase(TestCase):

    def setUp(self):
        self.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        self.variant = slowly_create_test_variant("3", 128198980, 'A', 'T', self.grch37)
        self.ref_variant = slowly_create_test_variant("3", 128198980, 'A', 'A', self.grch37)

    def test_variant_string(self):
        variant_string = "3:128198980 A>T"
        v = Variant.get_from_string(variant_string, self.grch37)
        self.assertEqual(self.variant, v)

    def test_variant_coordinate(self):
        vc = self.variant.coordinate
        v = Variant.get_from_variant_coordinate(vc, self.grch37)
        self.assertEqual(self.variant, v)

    def test_ref_variant_coordinate_explicit(self):
        vc_explicit = self.ref_variant.coordinate.as_external_explicit(self.grch37)
        self.assertEqual(vc_explicit.ref, vc_explicit.alt)

    def test_ref_variant_coordinate_internal(self):
        vc_symbolic = self.ref_variant.coordinate.as_internal_symbolic()
        self.assertEqual(vc_symbolic.alt, Variant.REFERENCE_ALT)
