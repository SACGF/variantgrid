"""GeneSymbol.resolve_name and GeneSymbol.resolve_names: which current symbol a name as written stands for."""
from django.test import TestCase

from genes.models import HGNC, GeneSymbol, GeneSymbolAlias, HGNCImport
from genes.models.models_gene import AmbiguousGeneSymbolAlias
from genes.models_enums import GeneSymbolAliasSource, HGNCStatus

# name as written -> (gene symbol, HGNC pk), or the error
EXPECTED = {
    "sept2": ("SEPTIN2", 7729),  # A rename outranks another gene's alias of the same name
    "ACPP": ("ACP3", 125),  # Previous symbol resolves though a GeneSymbol row shadows it
    "SEPTX": AmbiguousGeneSymbolAlias,  # An alias of several genes
    "AIP": ("AIP", 358),  # An approved symbol wins over another gene's informal alias
    "POPDC1": ("POPDC1", 1152),  # An approved symbol wins over a leftover redirect to its old name
    "bves": ("POPDC1", 1152),
    " SEPTIN6 ": ("SEPTIN6", 15848),
    "kal1": ("ANOS1", 3730),  # Alias stored in mixed case
    "RETIRED_ONLY": ("GONE", None),  # Alias of a symbol HGNC no longer has
    "NOT_A_GENE_SYMBOL": (None, None),
}


class TestGeneSymbolResolveName(TestCase):

    @classmethod
    def setUpTestData(cls):
        """ As the HGNC import stores them: each previous and alias symbol is also a GeneSymbolAlias """
        hgnc_import = HGNCImport.objects.create()

        def _hgnc(pk, symbol, previous_symbols="", alias_symbols=""):
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import, status=HGNCStatus.APPROVED,
                                approved_name=f"{symbol} approved name",
                                previous_symbols=previous_symbols, alias_symbols=alias_symbols)

        def _alias(alias, symbol, source=GeneSymbolAliasSource.HGNC):
            GeneSymbol.objects.get_or_create(symbol=symbol)
            GeneSymbolAlias.objects.create(alias=alias, gene_symbol_id=symbol, source=source)

        # SEPT2 is SEPTIN2's previous symbol and, on an unrelated gene, one of SEPTIN6's aliases
        _hgnc(7729, "SEPTIN2", previous_symbols="SEPT2")
        _hgnc(15848, "SEPTIN6", previous_symbols="SEPT6", alias_symbols="SEPT2")
        _alias("SEPT2", "SEPTIN2")
        _alias("SEPT2", "SEPTIN6")
        _alias("SEPTX", "SEPTIN2", GeneSymbolAliasSource.NCBI)
        _alias("SEPTX", "SEPTIN6", GeneSymbolAliasSource.NCBI)
        # Ensembl still has a GeneSymbol for ACP3's previous symbol
        GeneSymbol.objects.create(symbol="ACPP")
        _hgnc(125, "ACP3", previous_symbols="ACPP")
        _alias("ACPP", "ACP3")
        # AURKAIP1 lists the approved symbol AIP as an informal synonym
        _hgnc(358, "AIP")
        _hgnc(24114, "AURKAIP1", alias_symbols="AIP")
        _alias("AIP", "AURKAIP1")
        # A leftover row from when BVES was approved still points POPDC1 -> BVES
        _hgnc(1152, "POPDC1", previous_symbols="BVES")
        _alias("BVES", "POPDC1")
        _alias("POPDC1", "BVES")
        _hgnc(3730, "ANOS1", previous_symbols="KAL1")
        _alias("Kal1", "ANOS1")
        _alias("RETIRED_ONLY", "GONE")

    def _check(self, name, gene_symbol_id, hgnc):
        self.assertEqual(EXPECTED[name], (gene_symbol_id, hgnc.pk if hgnc else None), name)

    def test_resolve_name(self):
        for name, expected in EXPECTED.items():
            if expected is AmbiguousGeneSymbolAlias:
                with self.assertRaises(AmbiguousGeneSymbolAlias) as cm:
                    GeneSymbol.resolve_name(name)
                self.assertEqual(["SEPTIN2", "SEPTIN6"], cm.exception.gene_symbol_ids)
            else:
                self._check(name, *GeneSymbol.resolve_name(name))

    def test_resolve_names_together(self):
        with self.assertNumQueries(3):
            resolutions = GeneSymbol.resolve_names(EXPECTED)
        for name, (resolution, ambiguous) in resolutions.items():
            if EXPECTED[name] is AmbiguousGeneSymbolAlias:
                self.assertEqual(["SEPTIN2", "SEPTIN6"], ambiguous.gene_symbol_ids)
                self.assertIsNone(resolution.gene_symbol_id)
            else:
                self.assertIsNone(ambiguous)
                self._check(name, resolution.gene_symbol_id, resolution.hgnc)

    def test_resolve_names_records_alias_followed(self):
        resolutions = GeneSymbol.resolve_names(["ACPP", "AIP", "RETIRED_ONLY"])
        acpp_alias = GeneSymbolAlias.objects.get(alias="ACPP", gene_symbol_id="ACP3")
        retired_alias = GeneSymbolAlias.objects.get(alias="RETIRED_ONLY")
        self.assertEqual(acpp_alias.pk, resolutions["ACPP"][0].gene_symbol_alias_id)
        self.assertIsNone(resolutions["AIP"][0].gene_symbol_alias_id)
        self.assertEqual(retired_alias.pk, resolutions["RETIRED_ONLY"][0].gene_symbol_alias_id)
