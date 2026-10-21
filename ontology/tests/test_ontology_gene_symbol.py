"""OntologyTerm.get_gene_symbol: the HGNC term a gene name as written stands for."""
from django.test import TestCase
from django.utils import timezone

from genes.enums import GeneSymbolAliasSource, HGNCStatus
from genes.models import HGNC, GeneSymbol, GeneSymbolAlias, HGNCImport
from ontology.models import OntologyImport, OntologyService, OntologyTerm


class TestOntologyTermGetGeneSymbol(TestCase):

    @classmethod
    def setUpTestData(cls):
        """ SEPT2 is SEPTIN2's previous symbol and one of SEPTIN6's aliases, so both terms list it """
        hgnc_import = HGNCImport.objects.create()
        o_import = OntologyImport.objects.create(import_source=OntologyService.HGNC, filename="test",
                                                 context="test", hash="N/A", processor_version=1,
                                                 completed=True, processed_date=timezone.now())
        for pk, symbol, previous, alias in [(15848, "SEPTIN6", "SEPT6", "SEPT2"),
                                            (7729, "SEPTIN2", "SEPT2", "")]:
            GeneSymbol.objects.create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name",
                                previous_symbols=previous, alias_symbols=alias)
            GeneSymbolAlias.objects.create(alias="SEPT2", gene_symbol_id=symbol, source=GeneSymbolAliasSource.HGNC)
            OntologyTerm.objects.create(id=f"HGNC:{pk}", ontology_service=OntologyService.HGNC, index=pk,
                                        name=symbol, aliases=["SEPT2"], from_import=o_import)

    def test_a_name_on_several_terms_goes_to_the_gene_it_was_renamed_from(self):
        self.assertEqual("HGNC:7729", OntologyTerm.get_gene_symbol("SEPT2").id)
