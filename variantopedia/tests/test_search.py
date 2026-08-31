from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from annotation.fake_annotation import get_fake_annotation_version
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.models import Gene, GeneSymbol, Transcript, TranscriptVersion
from snpdb.models import ClinGenAllele, GenomeBuild, Quad, Trio, Variant
from snpdb.search import search_data
from snpdb.tests.utils.fake_cohort_data import create_fake_quad, create_fake_trio
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

VARIANT_CHROM = "1"
VARIANT_POSITION = 169519049
CLINGEN_ALLELE_ID = 285410130


@override_settings(PREFER_ALLELE_LINKS=False)
class TestSearch(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser')[0]

        get_fake_annotation_version(GenomeBuild.grch37())
        get_fake_annotation_version(GenomeBuild.grch38())

        create_fake_transcript_version(GenomeBuild.grch38())
        # Searches only return a Variant that's actually in the DB - without one every
        # coordinate/HGVS search below would pass over an empty result set
        cls.variant = slowly_create_test_variant(VARIANT_CHROM, VARIANT_POSITION, "T", "C",
                                                 GenomeBuild.grch37())

    def _verify_all_of_type(self, search_results, search_type, ignore_errors=False):
        checked = 0
        for sr in search_results.results:
            if ignore_errors:
                if sr.preview.is_error:
                    continue
            self.assertEqual(search_type, sr.search_type, f"Search result {sr} is of type {search_type}")
            checked += 1
        self.assertTrue(checked, f"No {search_type} results to check")

    def _assert_found_of_type(self, search_results, obj, search_type):
        matches = [sr for sr in search_results.results if sr.search_type == search_type]
        self.assertIn(obj, [sr.preview.obj for sr in matches],
                      f"{obj} not found in {search_type} search results")

    def test_search_hgvs(self):
        HGVS_NAMES = [
            "ENST00000300305.7(RUNX1):c.352-1G>A",
            "NC_000007.13:g.117199563G>T",
            "NM_001754.5:557T>A",  # This is a sloppy HGVS (missing c.) that needs to be cleaned
        ]
        for hgvs_name in HGVS_NAMES:
            search_results = search_data(self.user, hgvs_name, False)
            # Need to ignore errors as 37 HGVS uses fake clingen
            self._verify_all_of_type(search_results, Variant.preview_category(), ignore_errors=True)

    def test_search_locus(self):
        """ A locus (no ref/alt) resolves to the variants at that position """
        LOCI = [
            f"chr{VARIANT_CHROM}:{VARIANT_POSITION}",
            f"{VARIANT_CHROM}:{VARIANT_POSITION}",
        ]
        for locus_str in LOCI:
            search_results = search_data(self.user, locus_str, False)
            self._assert_found_of_type(search_results, self.variant, Variant.preview_category())

    def test_search_variant(self):
        VARIANTS = [
            f"{VARIANT_CHROM}:{VARIANT_POSITION} T>C",
            f"chr{VARIANT_CHROM}:{VARIANT_POSITION} T>C",
            f"{VARIANT_CHROM}-{VARIANT_POSITION}-T-C",  # gnomAD style
        ]
        for variant_str in VARIANTS:
            search_results = search_data(self.user, variant_str, False)
            # Returns variant if settings.PREFER_ALLELE_LINKS=False
            self._assert_found_of_type(search_results, self.variant, Variant.preview_category())

    def test_clingen_allele(self):
        """ A ClinGen allele resolves via its genomic HGVS to the variant in this build """
        api_response = {
            "genomicAlleles": [
                {"referenceGenome": "GRCh37", "chromosome": VARIANT_CHROM,
                 "hgvs": [f"NC_000001.10:g.{VARIANT_POSITION}T>C"]},
            ],
        }
        ClinGenAllele.objects.get_or_create(id=CLINGEN_ALLELE_ID, defaults={"api_response": api_response})
        search_results = search_data(self.user, f"CA{CLINGEN_ALLELE_ID}", False)
        self._assert_found_of_type(search_results, self.variant, Variant.preview_category())

    def test_gene_symbol(self):
        GENE_SYMBOL = [
            "runx1",
            "RUNX1",
        ]
        for gene_symbol in GENE_SYMBOL:
            search_results = search_data(self.user, gene_symbol, False)
            self._verify_all_of_type(search_results, GeneSymbol.preview_category())

    def test_gene(self):
        GENE = [
            "ENSG00000159216",
        ]
        for gene in GENE:
            search_results = search_data(self.user, gene, False)
            self._verify_all_of_type(search_results, Gene.preview_category())

    def test_transcript(self):
        TRANSCRIPTS = [
            "ENST00000300305",
        ]
        for transcript_id in TRANSCRIPTS:
            search_results = search_data(self.user, transcript_id, False)
            self._verify_all_of_type(search_results, Transcript.preview_category())

    def test_transcript_version(self):
        TRANSCRIPT_VERSION = [
            "ENST00000300305.7",
        ]
        for transcript_version in TRANSCRIPT_VERSION:
            search_results = search_data(self.user, transcript_version, False)
            self._verify_all_of_type(search_results, TranscriptVersion.preview_category())

    def test_search_trio(self):
        trio = create_fake_trio(self.user, GenomeBuild.grch38())
        search_results = search_data(self.user, trio.name, False)
        self._assert_found_of_type(search_results, trio, Trio.preview_category())

    def test_search_quad(self):
        quad = create_fake_quad(self.user, GenomeBuild.grch38())
        search_results = search_data(self.user, quad.name, False)
        self._assert_found_of_type(search_results, quad, Quad.preview_category())
