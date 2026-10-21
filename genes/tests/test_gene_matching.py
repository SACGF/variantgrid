from django.contrib.auth.models import User
from django.core.management import call_command
from django.test.testcases import TestCase

from genes.enums import AnnotationConsortium, GeneSymbolAliasSource, HGNCStatus
from genes.gene_matching import (
    MAX_GENE_SYMBOL_LENGTH,
    GeneSymbolMatcher,
    ReleaseGeneMatcher,
    tokenize_gene_symbols,
)
from genes.models import (
    HGNC,
    Gene,
    GeneAnnotationImport,
    GeneAnnotationRelease,
    GeneList,
    GeneSymbol,
    GeneSymbolAlias,
    GeneVersion,
    HGNCImport,
    ReleaseGeneSymbol,
    ReleaseGeneSymbolGene,
    ReleaseGeneVersion,
)
from snpdb.models import GenomeBuild


class TestGeneSymbolMatcherGeneList(TestCase):
    def test_alias_recorded_and_ambiguous_alias_warned(self):
        for symbol in ["SEPTIN2", "SEPTIN6"]:
            GeneSymbol.objects.create(symbol=symbol)
            GeneSymbolAlias.objects.create(alias="SEPTX", gene_symbol_id=symbol, source=GeneSymbolAliasSource.NCBI)
        kal1 = GeneSymbolAlias.objects.create(alias="KAL1", gene_symbol=GeneSymbol.objects.create(symbol="ANOS1"),
                                              source=GeneSymbolAliasSource.HGNC)
        user = User.objects.get_or_create(username='test_user')[0]
        gene_list = GeneList.objects.create(name="test_ambiguous", user=user)

        GeneSymbolMatcher().create_gene_list_gene_symbols_from_text(gene_list, "SEPTX KAL1")

        glgs_by_name = {glgs.original_name: glgs for glgs in gene_list.genelistgenesymbol_set.all()}
        self.assertIsNone(glgs_by_name["SEPTX"].gene_symbol_id)
        self.assertEqual(("ANOS1", kal1.pk), (glgs_by_name["KAL1"].gene_symbol_id,
                                              glgs_by_name["KAL1"].gene_symbol_alias_id))
        gene_list.refresh_from_db()
        self.assertIn("'SEPTX' is an ambiguous alias for SEPTIN2, SEPTIN6", gene_list.error_message)


class TestSetGeneListFromText(TestCase):
    def test_diff_keeps_alias_and_removes_unmatched(self):
        for symbol in ["ANOS1", "RUNX1", "GATA2"]:
            GeneSymbol.objects.create(symbol=symbol)
        GeneSymbolAlias.objects.create(alias="KAL1", gene_symbol_id="ANOS1", source=GeneSymbolAliasSource.HGNC)
        user = User.objects.get_or_create(username='test_user')[0]
        gene_list = GeneList.objects.create(name="test_set_from_text", user=user)
        matcher = GeneSymbolMatcher()
        matcher.create_gene_list_gene_symbols(gene_list, ["KAL1", "RUNX1", "NOTAGENE"], modification_info="original")

        # Text box shows the alias row as its symbol, and the unmatched row by its name
        added, removed = matcher.set_gene_list_gene_symbols_from_text(gene_list, "anos1, GATA2", modification_info="edit")

        self.assertEqual((["GATA2"], ["NOTAGENE", "RUNX1"]), (added, removed))
        rows = dict(gene_list.genelistgenesymbol_set.values_list("original_name", "modification_info"))
        self.assertEqual({"KAL1": "original", "GATA2": "edit"}, rows)

    def test_warning_follows_the_text(self):
        """ The ambiguous alias warning goes once the name is replaced, and stays while the name is kept """
        for symbol in ["SEPTIN2", "SEPTIN6", "RUNX1"]:
            GeneSymbol.objects.create(symbol=symbol)
            GeneSymbolAlias.objects.create(alias="SEPTX", gene_symbol_id=symbol, source=GeneSymbolAliasSource.NCBI)
        user = User.objects.get_or_create(username='test_user')[0]
        gene_list = GeneList.objects.create(name="test_warning", user=user)
        matcher = GeneSymbolMatcher()
        matcher.create_gene_list_gene_symbols_from_text(gene_list, "SEPTX")
        gene_list.refresh_from_db()
        self.assertIn("'SEPTX' is an ambiguous alias", gene_list.error_message)

        matcher.set_gene_list_gene_symbols_from_text(gene_list, "SEPTX RUNX1")
        gene_list.refresh_from_db()
        self.assertIn("'SEPTX' is an ambiguous alias", gene_list.error_message)

        matcher.set_gene_list_gene_symbols_from_text(gene_list, "SEPTIN2 RUNX1")
        gene_list.refresh_from_db()
        self.assertFalse(gene_list.error_message)


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


class ReleaseGeneMatcherTestCase(TestCase):
    """ Base fixture: a release with PDCD2/ANOS1/OLDNAME, and HGNC aliases where "RP8" is shared by
        MT-TS2 (absent from the release) and PDCD2. """

    @classmethod
    def setUpTestData(cls):
        genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        cls.gene_annotation_import = GeneAnnotationImport.objects.create(
            url="fake", genome_build=genome_build, annotation_consortium=AnnotationConsortium.REFSEQ)
        cls.release = GeneAnnotationRelease.objects.create(
            version="test_1669", annotation_consortium=AnnotationConsortium.REFSEQ,
            genome_build=genome_build, gene_annotation_import=cls.gene_annotation_import)

        def _release_gene(symbol, gene_id):
            gene_symbol = GeneSymbol.objects.get_or_create(symbol=symbol)[0]
            gene = Gene.objects.create(identifier=gene_id, annotation_consortium=AnnotationConsortium.REFSEQ)
            gene_version = GeneVersion.objects.create(gene=gene, gene_symbol=gene_symbol, version=1,
                                                     genome_build=genome_build,
                                                     import_source=cls.gene_annotation_import)
            ReleaseGeneVersion.objects.create(release=cls.release, gene_version=gene_version)
            return gene

        _release_gene("PDCD2", "5134")
        _release_gene("ANOS1", "3730")
        _release_gene("OLDNAME", "111")  # Release still carries this gene under its previous symbol

        def _alias(alias, symbol):
            GeneSymbol.objects.get_or_create(symbol=alias)
            GeneSymbol.objects.get_or_create(symbol=symbol)
            GeneSymbolAlias.objects.create(alias=alias, gene_symbol_id=symbol,
                                           source=GeneSymbolAliasSource.HGNC)

        # HGNC has RP8 as an alias for both MT-TS2 (not in release) and PDCD2 (in release)
        _alias("RP8", "MT-TS2")
        _alias("RP8", "PDCD2")
        _alias("KAL1", "ANOS1")
        _alias("OLDNAME", "NEWNAME")


class TestReleaseGeneMatcher(ReleaseGeneMatcherTestCase):
    """ A symbol reaches a release's genes directly, or through exactly one GeneSymbolAlias hop (either
        direction). Chaining hops let an alias string shared by two unrelated genes bridge them (#1669). """

    def _match(self, symbol):
        gm = ReleaseGeneMatcher(self.release)
        return gm._get_gene_id_and_match_info_for_symbol([symbol])[symbol]

    def test_direct_symbol_wins(self):
        self.assertEqual([("5134", None)], self._match("PDCD2"))

    def test_forward_alias_hop(self):
        self.assertEqual([("3730", "KAL1 is an alias for ANOS1 (HGNC)")], self._match("KAL1"))

    def test_backward_alias_hop(self):
        self.assertEqual([("111", "OLDNAME is an alias for NEWNAME (HGNC)")], self._match("NEWNAME"))

    def test_shared_alias_does_not_bridge(self):
        self.assertEqual([], self._match("MT-TS2"))

    def test_alias_matches_only_release_target(self):
        self.assertEqual([("5134", "RP8 is an alias for PDCD2 (HGNC)")], self._match("RP8"))


class TestFixRematchReleaseSymbolsToGenes(ReleaseGeneMatcherTestCase):
    """ The resync command has to delete matches the current rules no longer make - a rematch only inserts. """

    def setUp(self):
        release_gene_symbols = {}
        for symbol in ["PDCD2", "KAL1", "MT-TS2"]:
            release_gene_symbols[symbol] = ReleaseGeneSymbol.objects.create(release=self.release,
                                                                            gene_symbol_id=symbol)
        # What the old multi-hop traversal wrote: MT-TS2 -> RP8 -> PDCD2
        ReleaseGeneSymbolGene.objects.create(release_gene_symbol=release_gene_symbols["MT-TS2"], gene_id="5134",
                                             match_info="RP8 is an alias for PDCD2 (HGNC)")
        ReleaseGeneSymbolGene.objects.create(release_gene_symbol=release_gene_symbols["PDCD2"], gene_id="5134")

    def _matched_genes(self, symbol):
        return set(ReleaseGeneSymbolGene.objects.filter(release_gene_symbol__release=self.release,
                                                        release_gene_symbol__gene_symbol_id=symbol)
                   .values_list("gene_id", flat=True))

    def test_resync_deletes_chained_inserts_missing_and_keeps_direct(self):
        call_command("fix_rematch_release_symbols_to_genes")
        self.assertEqual(set(), self._matched_genes("MT-TS2"))
        self.assertEqual({"5134"}, self._matched_genes("PDCD2"))
        self.assertEqual({"3730"}, self._matched_genes("KAL1"))

    def test_dry_run_writes_nothing(self):
        call_command("fix_rematch_release_symbols_to_genes", dry_run=True)
        self.assertEqual({"5134"}, self._matched_genes("MT-TS2"))
        self.assertEqual(set(), self._matched_genes("KAL1"))
