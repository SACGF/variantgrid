from django.contrib.auth.models import User
from django.test.testcases import TestCase

from genes.gene_matching import GeneSymbolMatcher, HGNCMatcher, MAX_GENE_SYMBOL_LENGTH, tokenize_gene_symbols
from genes.models import GeneList, GeneSymbol, GeneSymbolAlias, HGNC, HGNCImport
from genes.models_enums import GeneSymbolAliasSource, HGNCStatus


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


class TestHGNCMatcher(TestCase):
    """ HGNC alias rows accumulate across imports and are never pruned, and include HGNC's informal
        alias_symbol synonyms - which are often another gene's approved symbol. So a symbol that HGNC
        currently approves has to win over any alias redirect. """

    def setUp(self):
        hgnc_import = HGNCImport.objects.create()

        def _hgnc(pk, symbol, previous_symbols="", alias_symbols=""):
            GeneSymbol.objects.get_or_create(symbol=symbol)
            return HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                       status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name",
                                       previous_symbols=previous_symbols, alias_symbols=alias_symbols)

        def _alias(alias, symbol):
            GeneSymbol.objects.get_or_create(symbol=symbol)
            GeneSymbolAlias.objects.create(alias=alias, gene_symbol_id=symbol,
                                           source=GeneSymbolAliasSource.HGNC)

        # AIP is approved, but AURKAIP1 lists "AIP" as an informal synonym
        self.aip = _hgnc(358, "AIP")
        self.aurkaip1 = _hgnc(24114, "AURKAIP1", alias_symbols="AIP")
        _alias("AIP", "AURKAIP1")

        # POPDC1 is approved after being renamed from BVES, but a leftover row from when BVES was
        # the approved symbol still points POPDC1 -> BVES, which no longer exists
        self.popdc1 = _hgnc(1152, "POPDC1", previous_symbols="BVES")
        _alias("BVES", "POPDC1")
        _alias("POPDC1", "BVES")

    def test_approved_symbol_wins_over_informal_alias(self):
        self.assertEqual(self.aip, HGNCMatcher().match_hgnc("AIP"))

    def test_approved_symbol_wins_over_retired_redirect(self):
        self.assertEqual(self.popdc1, HGNCMatcher().match_hgnc("POPDC1"))

    def test_previous_symbol_still_resolves_via_alias(self):
        self.assertEqual(self.popdc1, HGNCMatcher().match_hgnc("BVES"))

    def test_alias_is_case_insensitive(self):
        self.assertEqual(self.popdc1, HGNCMatcher().match_hgnc("bves"))

    def test_unknown_symbol_returns_none(self):
        self.assertIsNone(HGNCMatcher().match_hgnc("NOT_A_GENE_SYMBOL"))

    def test_alias_target_that_no_longer_exists_returns_none(self):
        _ = GeneSymbol.objects.get_or_create(symbol="GONE")
        GeneSymbolAlias.objects.create(alias="RETIRED_ONLY", gene_symbol_id="GONE",
                                       source=GeneSymbolAliasSource.HGNC)
        self.assertIsNone(HGNCMatcher().match_hgnc("RETIRED_ONLY"))


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
