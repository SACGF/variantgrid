"""classification_export_user_strings_to_q: the gene names a user types to restrict an export."""
from django.test import TestCase

from classification.views.exports.classification_export_filter import (
    classification_export_user_strings_to_q,
)
from genes.enums import GeneSymbolAliasSource
from genes.models import GeneSymbol, GeneSymbolAlias
from snpdb.models import GenomeBuild


class TestClassificationExportUserStrings(TestCase):

    @classmethod
    def setUpTestData(cls):
        for symbol in ["RUNX1", "BRCA1"]:
            GeneSymbol.objects.create(symbol=symbol)
            GeneSymbolAlias.objects.create(alias="SHARED", gene_symbol_id=symbol, source=GeneSymbolAliasSource.NCBI)
        GeneSymbolAlias.objects.create(alias="AML1", gene_symbol_id="RUNX1", source=GeneSymbolAliasSource.HGNC)

    def test_an_alias_filters_on_its_gene(self):
        q = classification_export_user_strings_to_q(["AML1"], GenomeBuild.grch37())
        self.assertIn("RUNX1", str(q))

    def test_an_alias_of_several_genes_is_an_error_naming_them(self):
        with self.assertRaisesMessage(ValueError, "'SHARED' is an ambiguous alias for BRCA1, RUNX1"):
            classification_export_user_strings_to_q(["SHARED"], GenomeBuild.grch37())
