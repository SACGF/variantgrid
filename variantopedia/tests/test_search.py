from unittest import TestCase

from django.contrib.auth.models import User

from variantopedia.search import search_data


class TestVCFProcessors(TestCase):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.user = User.objects.get_or_create(username='testuser')[0]

    def test_search_hgvs(self):
        HGVS_NAMES = [
            "NM_001080463.1:c.5972T>A",
            "NM_000492.3(CFTR):c.1438G>T",
            "NC_000007.13:g.117199563G>T",
        ]
        for hgvs_name in HGVS_NAMES:
            _search_results = search_data(self.user, hgvs_name, False)

    def test_search_gene_hgvs(self):
        HGVS_NAMES = [
            "CFTR:c.350G>A",
        ]
        for hgvs_name in HGVS_NAMES:
            _search_results = search_data(self.user, hgvs_name, False)

    def test_search_locus(self):
        LOCI = [
            "chr1:169519049",
        ]
        for locus_str in LOCI:
            _search_results = search_data(self.user, locus_str, False)

    def test_search_variant(self):
        VARIANTS = [
            "1:169519049 T>C",
            "chr1:169519049 T>C",
            "1-169519049-T-C",  # gnomAD style
        ]
        for variant_str in VARIANTS:
            _search_results = search_data(self.user, variant_str, False)

    def test_search_dbsnp(self):
        DBSNP = [
            "rs6025",
        ]
        for dbsnp_str in DBSNP:
            _search_results = search_data(self.user, dbsnp_str, False)

    def test_clingen_allele(self):
        CLINGEN_ALLELE = [
            "CA285410130",
        ]
        for ca_str in CLINGEN_ALLELE:
            _search_results = search_data(self.user, ca_str, False)

    def test_gene_symbol(self):
        GENE_SYMBOL = [
            "gata2",
            "GATA2",
        ]
        for gene_symbol in GENE_SYMBOL:
            _search_results = search_data(self.user, gene_symbol, False)

    def test_gene(self):
        GENE = [
            "ENSG00000179348",
        ]
        for gene in GENE:
            _search_results = search_data(self.user, gene, False)

    def test_transcript(self):
        TRANSCRIPTS = [
            "NM_001080463",
        ]
        for transcript_id in TRANSCRIPTS:
            _search_results = search_data(self.user, transcript_id, False)

    def test_transcript_version(self):
        TRANSCRIPT_VERSION = [
            "NM_001080463.2",
        ]
        for transcript_version in TRANSCRIPT_VERSION:
            _search_results = search_data(self.user, transcript_version, False)



