from unittest.mock import patch

from django.test import TestCase

from genes.models import GeneSymbol, PanelAppServer, PanelAppPanel, PanelAppPanelLocalCacheGeneSymbol, HGNC, HGNCImport
from genes.models_enums import HGNCStatus
from genes.panel_app import get_panel_app_local_cache, get_panel_app_panel_as_gene_list_json


class TestPanelAppLocalCache(TestCase):
    def setUp(self):
        self.server = PanelAppServer.objects.create(name="Test PanelApp",
                                                     url="https://panelapp.example.com",
                                                     icon_css_class="")
        self.panel = PanelAppPanel.objects.create(server=self.server,
                                                  panel_id=1,
                                                  disease_group="",
                                                  disease_sub_group="",
                                                  name="Test Panel",
                                                  status="public",
                                                  current_version="0.1")
        hgnc_import = HGNCImport.objects.create()
        # PanelApp reports the pre-rename symbol for this gene, we hold the current approved one
        GeneSymbol.objects.get_or_create(symbol="CFAP276")
        self.hgnc = HGNC.objects.create(pk=32331, gene_symbol_id="CFAP276", hgnc_import=hgnc_import,
                                        status=HGNCStatus.APPROVED, approved_name="cilia and flagella associated 276",
                                        previous_symbols="C1orf194", alias_symbols="")

    def _api_json(self, gene_data_list):
        return {
            "id": self.panel.panel_id,
            "name": self.panel.name,
            "disease_group": "",
            "disease_sub_group": "",
            "status": "public",
            "version": "0.2",
            "relevant_disorders": [],
            "genes": [{"gene_data": gene_data, "confidence_level": "3"} for gene_data in gene_data_list],
        }

    def _cache_panel(self, mock_api_json, gene_data_list):
        mock_api_json.return_value = self._api_json(gene_data_list)
        return get_panel_app_local_cache(self.panel)

    @patch("genes.panel_app._get_panel_app_panel_api_json")
    def test_resolves_via_hgnc_id_not_symbol(self, mock_api_json):
        """ PanelApp's gene symbol comes from a dated HGNC snapshot, so a stale symbol paired with a
            current HGNC ID must resolve to our approved symbol. """
        pap_lc = self._cache_panel(mock_api_json, [{"gene_symbol": "C1orf194", "hgnc_id": "HGNC:32331"}])

        record = PanelAppPanelLocalCacheGeneSymbol.objects.get(panel_app_local_cache=pap_lc)
        self.assertEqual(self.hgnc, record.hgnc)
        self.assertEqual("C1orf194", record.gene_symbol_reported)
        self.assertEqual("CFAP276", record.gene_symbol_str)

    @patch("genes.panel_app._get_panel_app_panel_api_json")
    def test_unknown_hgnc_id_falls_back_to_reported_symbol(self, mock_api_json):
        pap_lc = self._cache_panel(mock_api_json, [{"gene_symbol": "NOT_IN_HGNC", "hgnc_id": "HGNC:99999999"}])

        record = PanelAppPanelLocalCacheGeneSymbol.objects.get(panel_app_local_cache=pap_lc)
        self.assertIsNone(record.hgnc)
        self.assertEqual("NOT_IN_HGNC", record.gene_symbol_str)

    @patch("genes.panel_app._get_panel_app_panel_api_json")
    def test_missing_hgnc_id_falls_back_to_reported_symbol(self, mock_api_json):
        pap_lc = self._cache_panel(mock_api_json, [{"gene_symbol": "BRCA1"}])

        record = PanelAppPanelLocalCacheGeneSymbol.objects.get(panel_app_local_cache=pap_lc)
        self.assertIsNone(record.hgnc)
        self.assertEqual("BRCA1", record.gene_symbol_str)

    @patch("genes.panel_app._get_panel_app_panel_api_json")
    def test_unparsable_hgnc_id_falls_back_to_reported_symbol(self, mock_api_json):
        pap_lc = self._cache_panel(mock_api_json, [{"gene_symbol": "BRCA1", "hgnc_id": "HGNC:not-a-number"}])

        record = PanelAppPanelLocalCacheGeneSymbol.objects.get(panel_app_local_cache=pap_lc)
        self.assertIsNone(record.hgnc)
        self.assertEqual("BRCA1", record.gene_symbol_str)

    @patch("genes.panel_app._get_panel_app_panel_api_json")
    def test_caching_does_not_create_gene_symbols_for_panel_app_names(self, mock_api_json):
        """ PanelApp names we don't recognise no longer mint GeneSymbol rows - a symbol from their
            older snapshot isn't one of ours. """
        self._cache_panel(mock_api_json, [{"gene_symbol": "C1orf194", "hgnc_id": "HGNC:32331"}])

        self.assertFalse(GeneSymbol.objects.filter(symbol="C1orf194").exists())

    @patch("genes.panel_app._get_panel_app_panel_api_json")
    def test_panel_preview_uses_our_approved_symbol(self, mock_api_json):
        """ The panel preview matched on PanelApp's symbol, so a gene they still call by its old name
            was dropped as unmatched. Resolving via HGNC ID keys the evidence by our symbol. """
        mock_api_json.return_value = self._api_json([{"gene_symbol": "C1orf194", "hgnc_id": "HGNC:32331"}])

        data = get_panel_app_panel_as_gene_list_json(self.panel.pk)

        self.assertEqual(["CFAP276"], list(data["gene_evidence"]))

    @patch("genes.panel_app._get_panel_app_panel_api_json")
    def test_gene_list_uses_our_approved_symbol(self, mock_api_json):
        pap_lc = self._cache_panel(mock_api_json, [{"gene_symbol": "C1orf194", "hgnc_id": "HGNC:32331"}])

        gene_list = pap_lc.get_gene_list(panel_app_confidence=3)
        symbols = set(gene_list.genelistgenesymbol_set.values_list("gene_symbol_id", flat=True))
        self.assertEqual({"CFAP276"}, symbols)
