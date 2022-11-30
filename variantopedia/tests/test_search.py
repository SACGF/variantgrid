from django.contrib.auth.models import User
from django.test import TestCase

from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from snpdb.models import GenomeBuild, ClinGenAllele
from variantopedia.models import SearchTypes
from variantopedia.search import search_data


class TestSearch(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser')[0]
        create_fake_transcript_version(GenomeBuild.grch38())

    def _verify_all_of_type(self, search_results, search_type):
        for sr in search_results.results:
            self.assertEqual(search_type, sr.search_type, f"Search result {sr} is of type {search_type}")

    def test_search_hgvs(self):
        HGVS_NAMES = [
            "ENST00000300305.7(RUNX1):c.352-1G>A",
            "NC_000007.13:g.117199563G>T",
            "NM_001754.5:557T>A",  # This is a sloppy HGVS (missing c.) that needs to be cleaned
        ]
        for hgvs_name in HGVS_NAMES:
            search_results = search_data(self.user, hgvs_name, False)
            self._verify_all_of_type(search_results, SearchTypes.HGVS)

    def test_search_gene_hgvs(self):
        HGVS_NAMES = [
            # "CFTR:c.50G>A",  # This is really slow...
        ]
        for hgvs_name in HGVS_NAMES:
            _search_results = search_data(self.user, hgvs_name, False)

    def test_search_locus(self):
        LOCI = [
            "chr1:169519049",
        ]
        for locus_str in LOCI:
            search_results = search_data(self.user, locus_str, False)
            self._verify_all_of_type(search_results, SearchTypes.LOCUS)

    def test_search_variant(self):
        VARIANTS = [
            "1:169519049 T>C",
            "chr1:169519049 T>C",
            "chrM:9049 T>C",
            "MT:9049 T>C",
            "1-169519049-T-C",  # gnomAD style
        ]
        for variant_str in VARIANTS:
            search_results = search_data(self.user, variant_str, False)
            self._verify_all_of_type(search_results, SearchTypes.VARIANT)

    def test_search_dbsnp(self):
        DBSNP = [
            # "rs6025",
        ]
        for dbsnp_str in DBSNP:
            _search_results = search_data(self.user, dbsnp_str, False)

    def test_clingen_allele(self):
        CLINGEN_ALLELE = [
            "CA285410130",
        ]
        ClinGenAllele.objects.get_or_create(id=285410130, defaults={"api_response": {}})
        for ca_str in CLINGEN_ALLELE:
            search_results = search_data(self.user, ca_str, False)
            self._verify_all_of_type(search_results, SearchTypes.CLINGEN_ALLELE_ID)

    def test_gene_symbol(self):
        GENE_SYMBOL = [
            "runx1",
            "RUNX1",
        ]
        for gene_symbol in GENE_SYMBOL:
            search_results = search_data(self.user, gene_symbol, False)
            self._verify_all_of_type(search_results, SearchTypes.GENE_SYMBOL)

    def test_gene(self):
        GENE = [
            "ENSG00000159216",
        ]
        for gene in GENE:
            search_results = search_data(self.user, gene, False)
            self._verify_all_of_type(search_results, SearchTypes.GENE)

    def test_transcript(self):
        TRANSCRIPTS = [
            "ENST00000300305",
        ]
        for transcript_id in TRANSCRIPTS:
            search_results = search_data(self.user, transcript_id, False)
            self._verify_all_of_type(search_results, SearchTypes.TRANSCRIPT)

    def test_transcript_version(self):
        TRANSCRIPT_VERSION = [
            "ENST00000300305.7",
        ]
        for transcript_version in TRANSCRIPT_VERSION:
            search_results = search_data(self.user, transcript_version, False)
            self._verify_all_of_type(search_results, SearchTypes.TRANSCRIPT)
