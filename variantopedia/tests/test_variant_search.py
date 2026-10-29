from types import SimpleNamespace

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from annotation.cosmic import CosmicAPI
from annotation.fake_data import get_fake_annotation_version
from genes.hgvs import HGVSMatcher
from library.preview_request import PreviewData
from snpdb.models import GenomeBuild, Sample, Sequence, Variant, VariantCoordinate
from snpdb.models.models_variant import VARIANT_PATTERN
from snpdb.search import (
    SearchInput,
    SearchMessageOverall,
    SearchResponse,
    SearchResponsesCombined,
    SearchResult,
    SearchResultMatchStrength,
)
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele, slowly_create_test_variant
from variantopedia.signals.variant_search import (
    _alt_description,
    _alt_mismatch_message,
    search_allele_id,
    search_variant_id,
    yield_search_variant_match,
)


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


class SearchAutoOpenGenomeBuildTestCase(TestCase):
    """ A single result outside the preferred build auto-opens only if it can never be in that build (#2174) """

    @classmethod
    def setUpTestData(cls):
        cls.search_input = SearchInput(user=None, search_string="x",
                                       genome_build_preferred=GenomeBuild.get_name_or_alias("GRCh38"))
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

    def _single_preferred_result(self, search_type):
        preview = PreviewData(category=search_type.preview_category(), identifier="x", internal_url="/result",
                              genome_builds={self.grch37}, can_liftover=search_type.preview_can_liftover())
        result = SearchResult(preview=preview, match_strength=SearchResultMatchStrength.STRONG_MATCH)
        response = SearchResponse(search_input=self.search_input, search_type=search_type, matched_pattern=True,
                                  results=[result], total_count=1)
        return SearchResponsesCombined(self.search_input, [response]).single_preferred_result()

    def test_fixed_build_result_in_other_build_auto_opens(self):
        self.assertIsNotNone(self._single_preferred_result(Sample))

    def test_variant_in_other_build_does_not_auto_open(self):
        self.assertIsNone(self._single_preferred_result(Variant))


@override_settings(SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN=True)
class SearchIdRequiresClassificationTestCase(TestCase):
    """ Searching by internal variant/allele id applies the same visibility rule as searching by coordinate """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user(username="search_id_non_admin")
        cls.admin = User.objects.create_superuser(username="search_id_admin")
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.variant = slowly_create_test_variant("3", 5000, "A", "T", cls.grch37)
        cls.allele = create_mock_allele(cls.variant, cls.grch37)

    def _result_count(self, user, receiver, search_string) -> int:
        search_input = SearchInput(user=user, search_string=search_string, genome_build_preferred=self.grch37)
        return len(receiver.search(search_input).results)

    def test_unclassified_variant_hidden_from_non_admin(self):
        self.assertEqual(self._result_count(self.user, search_variant_id, f"v{self.variant.pk}"), 0)
        self.assertEqual(self._result_count(self.admin, search_variant_id, f"v{self.variant.pk}"), 1)

    def test_unclassified_allele_hidden_from_non_admin(self):
        self.assertEqual(self._result_count(self.user, search_allele_id, f"a{self.allele.pk}"), 0)
        self.assertEqual(self._result_count(self.admin, search_allele_id, f"a{self.allele.pk}"), 1)


class VariantSearchUnknownContigTestCase(TestCase):
    """ A coordinate on a contig the build doesn't have (#2529) - T2T has no MT """

    def test_contig_not_in_build_is_reported_not_raised(self):
        t2t = GenomeBuild.get_name_or_alias("T2T-CHM13v2.0")
        search_input = SimpleNamespace(genome_builds=[t2t], match=VARIANT_PATTERN.match("MT:16519 T>C"))
        [response] = list(yield_search_variant_match(search_input, VariantCoordinate.from_variant_match))
        self.assertIsInstance(response, SearchMessageOverall)
        self.assertIn("'MT'", response.message)
