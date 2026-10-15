"""GeneSymbol.resolve_name: which current symbol a name as written stands for."""
from django.test import TestCase

from genes.models import HGNC, GeneSymbol, GeneSymbolAlias, HGNCImport
from genes.models.models_gene import AmbiguousGeneSymbolAlias
from genes.models_enums import GeneSymbolAliasSource, HGNCStatus


class TestGeneSymbolResolveName(TestCase):

    @classmethod
    def setUpTestData(cls):
        """ As the HGNC import stores them: each previous and alias symbol is also a GeneSymbolAlias """
        hgnc_import = HGNCImport.objects.create()
        # SEPT2 is SEPTIN2's previous symbol and, on an unrelated gene, one of SEPTIN6's aliases
        for pk, symbol, previous, alias in [(7729, "SEPTIN2", "SEPT2", ""),
                                            (15848, "SEPTIN6", "SEPT6", "SEPT2")]:
            GeneSymbol.objects.create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name",
                                previous_symbols=previous, alias_symbols=alias)
            GeneSymbolAlias.objects.create(alias="SEPT2", gene_symbol_id=symbol, source=GeneSymbolAliasSource.HGNC)
        # Ensembl still has a GeneSymbol for ACP3's previous symbol
        GeneSymbol.objects.create(symbol="ACP3")
        GeneSymbol.objects.create(symbol="ACPP")
        HGNC.objects.create(pk=125, gene_symbol_id="ACP3", hgnc_import=hgnc_import, status=HGNCStatus.APPROVED,
                            approved_name="acid phosphatase 3", previous_symbols="ACPP")
        GeneSymbolAlias.objects.create(alias="ACPP", gene_symbol_id="ACP3", source=GeneSymbolAliasSource.HGNC)

    def test_a_rename_outranks_another_gene_alias_of_the_same_name(self):
        symbol, hgnc = GeneSymbol.resolve_name("sept2")
        self.assertEqual("SEPTIN2", symbol)
        self.assertEqual(7729, hgnc.pk)

    def test_previous_symbol_resolves_when_a_gene_symbol_row_shadows_it(self):
        self.assertEqual("ACP3", GeneSymbol.resolve_name("ACPP")[0])

    def test_an_alias_of_several_genes_raises_naming_them(self):
        for symbol in ["SEPTIN2", "SEPTIN6"]:
            GeneSymbolAlias.objects.create(alias="SEPTX", gene_symbol_id=symbol, source=GeneSymbolAliasSource.NCBI)
        with self.assertRaises(AmbiguousGeneSymbolAlias) as cm:
            GeneSymbol.resolve_name("SEPTX")
        self.assertEqual(["SEPTIN2", "SEPTIN6"], cm.exception.gene_symbol_ids)
