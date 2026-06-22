from unittest.mock import patch

from django.test import TestCase

from genes.models import GeneSymbol, PanelAppServer, PanelAppPanel, PanelAppPanelLocalCacheGeneSymbol
from genes.panel_app import get_panel_app_local_cache


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

    def _api_json(self, symbols):
        return {
            "id": self.panel.panel_id,
            "name": self.panel.name,
            "disease_group": "",
            "disease_sub_group": "",
            "status": "public",
            "version": "0.2",
            "relevant_disorders": [],
            "genes": [{"gene_data": {"gene_symbol": symbol}} for symbol in symbols],
        }

    @patch("genes.panel_app._get_panel_app_panel_api_json")
    def test_creates_new_gene_symbols(self, mock_api_json):
        """ Regression for #1567: new gene symbols must be bulk-created as GeneSymbol
            instances, not bare strings (AttributeError: 'str' object has no attribute 'pk'). """
        new_symbol = "BRCA1"
        self.assertFalse(GeneSymbol.objects.filter(symbol=new_symbol).exists())
        mock_api_json.return_value = self._api_json([new_symbol])

        pap_lc = get_panel_app_local_cache(self.panel)

        self.assertTrue(GeneSymbol.objects.filter(symbol=new_symbol).exists())
        cached_symbols = set(PanelAppPanelLocalCacheGeneSymbol.objects.filter(
            panel_app_local_cache=pap_lc).values_list("gene_symbol_id", flat=True))
        self.assertEqual({new_symbol}, cached_symbols)
