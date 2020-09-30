from django.test.testcases import TestCase

from genes.gene_matching import GeneSymbolMatcher
from genes.models import GeneSymbol


class TestGeneMatching(TestCase):
    # TODO: This isn't really testing anything anymore, maybe put in some aliases??
    GENE_SYMBOLS = ["AGRN", "BRCA1"]

    def setUp(self):
        for gene_symbol in TestGeneMatching.GENE_SYMBOLS:
            gene_symbol, _ = GeneSymbol.objects.get_or_create(symbol=gene_symbol)

    def test_gene_matcher(self):
        gene_matcher = GeneSymbolMatcher()

        for gene_name in TestGeneMatching.GENE_SYMBOLS:
            gene_symbol_id = gene_matcher.get_gene_symbol_id(gene_name)
            self.assertEqual(gene_name, gene_symbol_id)
            # print(f"{gene_symbol}/{transcript_id} => gene: {gene}, transcript: {transcript}")
