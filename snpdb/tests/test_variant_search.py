from django.test import TestCase

from annotation.cosmic import CosmicAPI
from genes.hgvs import HGVSMatcher
from snpdb.models import GenomeBuild, Sequence, VariantCoordinate
from snpdb.signals.variant_search import _alt_description, _alt_mismatch_message
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


class VariantSearchAltMismatchMessageTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.found_variant = slowly_create_test_variant("3", 128198980, "A", "T", cls.grch37)
        cls.searched = VariantCoordinate(chrom="3", position=128198980, ref="A", alt="G")

    def test_hgvs_search_names_both_as_g_hgvs(self):
        message = _alt_mismatch_message(self.searched, self.found_variant, HGVSMatcher.instance(self.grch37))
        self.assertEqual(message, 'No results for "NC_000003.11:g.128198980A>G", '
                                  'but found "NC_000003.11:g.128198980A>T" at the same position')

    def test_coordinate_search_describes_alts(self):
        message = _alt_mismatch_message(self.searched, self.found_variant, None)
        self.assertEqual(message, 'No results for alt "G", but found this using alt "T"')


class CosmicSearchGenomeBuildTestCase(TestCase):
    """ The clinicaltables COSMIC API only has GRCh37/38 data - other builds 500 out (issue #1783) """

    def test_t2t_is_not_supported(self):
        t2t = GenomeBuild.t2tv2()
        self.assertFalse(CosmicAPI.supports_genome_build(t2t))
        with self.assertRaises(ValueError):
            CosmicAPI("COSV53567516", t2t)

    def test_grch37_and_38_are_supported(self):
        for build_name in ["GRCh37", "GRCh38"]:
            self.assertTrue(CosmicAPI.supports_genome_build(GenomeBuild.get_name_or_alias(build_name)))
