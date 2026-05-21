from django.contrib.auth.models import User
from django.test.testcases import TestCase

from genes.gene_matching import GeneSymbolMatcher, MAX_GENE_SYMBOL_LENGTH, tokenize_gene_symbols
from genes.models import GeneList, GeneSymbol


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


class TestOversizedTokens(TestCase):
    def test_tokenize_splits_oversized(self):
        oversized = "X" * (MAX_GENE_SYMBOL_LENGTH + 50)
        result = tokenize_gene_symbols(f"BRCA1 RUNX1 {oversized}")
        self.assertEqual({"BRCA1", "RUNX1"}, result.valid)
        self.assertEqual([oversized], result.oversized)

    def test_tokenize_no_oversized(self):
        result = tokenize_gene_symbols("BRCA1, RUNX1")
        self.assertEqual({"BRCA1", "RUNX1"}, result.valid)
        self.assertEqual([], result.oversized)

    def test_create_gene_list_sets_warning_for_oversized(self):
        user = User.objects.get_or_create(username='test_user')[0]
        gene_list = GeneList.objects.create(name="test_oversized", user=user)
        oversized = "Y" * (MAX_GENE_SYMBOL_LENGTH + 1)
        gene_text = f"BRCA1 RUNX1 {oversized}"

        gene_matcher = GeneSymbolMatcher()
        gene_matcher.create_gene_list_gene_symbols_from_text(gene_list, gene_text)

        gene_list.refresh_from_db()
        self.assertIn("oversized", gene_list.error_message.lower())
        self.assertIn(str(MAX_GENE_SYMBOL_LENGTH), gene_list.error_message)
        # Valid symbols still got saved
        saved_names = set(gene_list.genelistgenesymbol_set.values_list("original_name", flat=True))
        self.assertEqual({"BRCA1", "RUNX1"}, saved_names)

    def test_create_gene_list_no_warning_when_all_valid(self):
        user = User.objects.get_or_create(username='test_user')[0]
        gene_list = GeneList.objects.create(name="test_clean", user=user)
        gene_matcher = GeneSymbolMatcher()
        gene_matcher.create_gene_list_gene_symbols_from_text(gene_list, "BRCA1 RUNX1")

        gene_list.refresh_from_db()
        self.assertFalse(gene_list.error_message)

    def test_create_gene_list_filters_oversized_from_prefiltered_list(self):
        """ REST API path passes a pre-built list; oversized entries should still be caught. """
        user = User.objects.get_or_create(username='test_user')[0]
        gene_list = GeneList.objects.create(name="test_rest", user=user)
        oversized = "Z" * (MAX_GENE_SYMBOL_LENGTH + 5)

        gene_matcher = GeneSymbolMatcher()
        gene_matcher.create_gene_list_gene_symbols(gene_list, ["BRCA1", oversized])

        gene_list.refresh_from_db()
        self.assertIn("oversized", gene_list.error_message.lower())
        saved_names = set(gene_list.genelistgenesymbol_set.values_list("original_name", flat=True))
        self.assertEqual({"BRCA1"}, saved_names)
