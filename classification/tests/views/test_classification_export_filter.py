"""classification_export_user_strings_to_q: the gene names a user types to restrict an export."""
from django.test import TestCase

from classification.views.exports.classification_export_filter import (
    classification_export_user_strings_to_q,
)
from genes.enums import GeneSymbolAliasSource, HGNCStatus
from genes.models import HGNC, GeneSymbol, GeneSymbolAlias, HGNCImport
from snpdb.models import GenomeBuild


class TestClassificationExportUserStrings(TestCase):

    @classmethod
    def setUpTestData(cls):
        for symbol in ["RUNX1", "BRCA1"]:
            GeneSymbol.objects.create(symbol=symbol)
            GeneSymbolAlias.objects.create(alias="SHARED", gene_symbol_id=symbol, source=GeneSymbolAliasSource.NCBI)
        GeneSymbolAlias.objects.create(alias="AML1", gene_symbol_id="RUNX1", source=GeneSymbolAliasSource.HGNC)
        for symbol in ["SLU7", "SRSF7"]:
            GeneSymbol.objects.create(symbol=symbol)
            GeneSymbolAlias.objects.create(alias="9G8", gene_symbol_id=symbol, source=GeneSymbolAliasSource.NCBI)
        # An older gene annotation import left a GeneSymbol for SEPTIN14's previous symbol
        for symbol in ["SEPT14", "SEPTIN14"]:
            GeneSymbol.objects.create(symbol=symbol)
        HGNC.objects.create(pk=33280, gene_symbol_id="SEPTIN14", hgnc_import=HGNCImport.objects.create(),
                            status=HGNCStatus.APPROVED, approved_name="septin 14", previous_symbols="SEPT14")
        GeneSymbolAlias.objects.create(alias="SEPT14", gene_symbol_id="SEPTIN14", source=GeneSymbolAliasSource.HGNC)

    def test_an_alias_filters_on_its_gene(self):
        q = classification_export_user_strings_to_q(["AML1"], GenomeBuild.grch37())
        self.assertIn("RUNX1", str(q))

    def test_an_alias_of_several_genes_is_an_error_naming_them(self):
        with self.assertRaisesMessage(ValueError, "'SHARED' is an ambiguous alias for BRCA1, RUNX1"):
            classification_export_user_strings_to_q(["SHARED"], GenomeBuild.grch37())

    def test_an_old_name_with_its_own_gene_symbol_filters_on_the_approved_symbol(self):
        q = classification_export_user_strings_to_q(["SEPT14"], GenomeBuild.grch37())
        self.assertIn("SEPTIN14", str(q))
        self.assertNotIn("'SEPT14'", str(q))

    def test_an_ambiguous_alias_starting_with_a_digit_is_an_error_naming_its_genes(self):
        with self.assertRaisesMessage(ValueError, "'9G8' is an ambiguous alias for SLU7, SRSF7"):
            classification_export_user_strings_to_q(["9G8"], GenomeBuild.grch37())
