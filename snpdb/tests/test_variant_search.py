from django.test import TestCase

from snpdb.models import GenomeBuild, Sequence, VariantCoordinate
from snpdb.signals.variant_search import _alt_description
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class VariantSearchAltDescriptionTestCase(TestCase):
    """ Regression tests for _alt_description (issue #3864).

        When a variant search has no direct hit it looks for "alt alts" - existing Variants at the same
        locus with a different alt - and describes each one. A Variant's alt is a Sequence (not a str), so
        describing it via Sequence.abbreviate(variant.alt) sliced the Sequence object and raised
        "'Sequence' object is not subscriptable" once the alt exceeded the abbreviation length. """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.long_alt = "G" + "A" * 30  # >20 chars, so abbreviate() slices - this is what crashed
        cls.long_alt_variant = slowly_create_test_variant("3", 128198995, "G", cls.long_alt, cls.grch37)
        cls.short_alt_variant = slowly_create_test_variant("3", 128198980, "A", "T", cls.grch37)

    def test_long_alt_on_variant_does_not_crash(self):
        """ Variant.alt is a Sequence - long alts must be abbreviated without subscripting the Sequence. """
        # alt is a Sequence, not a str - this is the condition that triggered the bug
        self.assertIsInstance(self.long_alt_variant.alt, Sequence)
        desc = _alt_description(self.long_alt_variant)
        self.assertIsInstance(desc, str)
        self.assertEqual(desc, Sequence.abbreviate(self.long_alt))

    def test_short_alt_on_variant(self):
        desc = _alt_description(self.short_alt_variant)
        self.assertEqual(desc, "T")

    def test_alt_on_variant_coordinate(self):
        """ VariantCoordinate.alt is a plain str - long alts are still abbreviated. """
        vc = VariantCoordinate(chrom="3", position=128198995, ref="G", alt=self.long_alt)
        self.assertEqual(_alt_description(vc), Sequence.abbreviate(self.long_alt))
