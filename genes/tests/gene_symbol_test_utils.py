"""Gene symbol / alias / HGNC rows for testing GeneSymbol.resolve_names and the pages that explain it."""
from genes.enums import GeneSymbolAliasSource, HGNCStatus
from genes.models import HGNC, GeneSymbol, GeneSymbolAlias, HGNCImport


def create_gene_symbol_alias_fixture():
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
