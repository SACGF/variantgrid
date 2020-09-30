from typing import Dict

import requests
from rest_framework.exceptions import NotFound

from genes.gene_matching import GeneSymbolMatcher
from genes.models import PanelAppPanelRelevantDisorders, PanelAppPanel, FakeGeneList
from genes.serializers import GeneListGeneSymbolSerializer


PANEL_APP_PREFIX = "panel-app-"
PANEL_APP_BASE_URL = "https://panelapp.genomicsengland.co.uk"
PANEL_APP_LIST_PANELS_URL = PANEL_APP_BASE_URL + "/WebServices/list_panels"
PANEL_APP_GET_PANEL_API_BASE_URL = PANEL_APP_BASE_URL + "/WebServices/get_panel/"
PANEL_APP_VIEW_PANEL_BASE_URL = PANEL_APP_BASE_URL + "/panels/"
PANEL_APP_SEARCH_BY_GENES_BASE_URL = PANEL_APP_BASE_URL + "/WebServices/search_genes/"


def get_panel_app_results_by_gene_symbol_json(gene_symbol):
    url = PANEL_APP_SEARCH_BY_GENES_BASE_URL + gene_symbol
    r = requests.get(url)
    data = r.json()
    return data.get("results")


def get_panel_app_panel_as_gene_list_json(gene_list_id, panel_app_panel_id):
    url = PANEL_APP_GET_PANEL_API_BASE_URL + panel_app_panel_id
    r = requests.get(url)
    json_data: Dict = r.json()
    # Panel App isn't very REST-ful - returns 200 for missing data but we'll return 404
    if isinstance(json_data, list):
        first_element = json_data[0]
        if first_element.endswith("not found."):
            detail = f"PanelApp couldn't find {panel_app_panel_id} ({url}) - returned {first_element}"
            raise NotFound(detail=detail)

    result = json_data["result"]
    genes = result["Genes"]
    name = result["SpecificDiseaseName"]
    absolute_url = PANEL_APP_VIEW_PANEL_BASE_URL + panel_app_panel_id

    gene_matcher = GeneSymbolMatcher()
    gene_list = FakeGeneList(name=name, user=None)
    gene_names_list = [gene_data["GeneSymbol"] for gene_data in genes]
    gene_list_gene_symbols = gene_matcher.create_gene_list_gene_symbols(gene_list, gene_names_list, save=False)
    genelistgenesymbol_set = GeneListGeneSymbolSerializer(gene_list_gene_symbols, many=True).data  # sorted(glg_set)

    # Need to build evidence from panel app but using our gene names
    panel_app_gene_evidence = {}
    for gene_data in genes:
        gene_symbol = gene_data["GeneSymbol"]
        panel_app_gene_evidence[gene_symbol] = gene_data

    gene_evidence = {}
    for glgs in gene_list_gene_symbols:
        if glgs.gene_symbol_id:
            evidence = panel_app_gene_evidence[glgs.original_name]
            gene_evidence[glgs.gene_symbol_id] = evidence
            # TODO: Handle unmatched symbols??

    data = {"pk": gene_list_id,
            "category": {"name": "PanelApp", "icon_css_class": "panel-app-icon"},
            "name": name,
            "import_status": "S",
            "genelistgenesymbol_set": genelistgenesymbol_set,
            "can_write": False,
            "absolute_url": absolute_url,
            "gene_evidence": gene_evidence}
    return data


def store_panel_app_panels_from_web(cached_web_resource):
    """ Used to build autocomplete list """

    r = requests.get(PANEL_APP_LIST_PANELS_URL)
    data = r.json()
    num_panels = 0

    for result in data["result"]:
        num_panels += 1
        pap = PanelAppPanel.objects.create(panel_id=result['Panel_Id'],
                                           cached_web_resource=cached_web_resource,
                                           disease_group=result['DiseaseGroup'],
                                           disease_sub_group=result['DiseaseSubGroup'],
                                           name=result['Name'],
                                           current_version=result['CurrentVersion'])

        relevant_disorders = result['Relevant_disorders']
        for rd in relevant_disorders:
            PanelAppPanelRelevantDisorders.objects.create(panel_app_panel=pap,
                                                          name=rd)

    cached_web_resource.description = f"{num_panels} panel app panels."
    cached_web_resource.save()
