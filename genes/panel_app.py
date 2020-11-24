import logging
from typing import Dict

import requests
from rest_framework.exceptions import NotFound

from genes.gene_matching import GeneSymbolMatcher
from genes.models import PanelAppPanelRelevantDisorders, PanelAppPanel, FakeGeneList
from genes.serializers import GeneListGeneSymbolSerializer


PANEL_APP_PREFIX = "panel-app-"
PANEL_APP_BASE_URL = "https://panelapp.genomicsengland.co.uk"
PANEL_APP_LIST_PANELS_PATH = "/api/v1/panels/"
PANEL_APP_GET_PANEL_API_BASE_PATH = "/api/v1/panels/"
PANEL_APP_SEARCH_BY_GENES_BASE_PATH = "/api/v1/genes/"


def get_panel_app_results_by_gene_symbol_json(gene_symbol):
    url = PANEL_APP_BASE_URL + PANEL_APP_SEARCH_BY_GENES_BASE_PATH + gene_symbol
    r = requests.get(url)
    data = r.json()
    return data.get("results")


def get_panel_app_panel_as_gene_list_json(gene_list_id, panel_app_panel_id):
    url = PANEL_APP_BASE_URL + PANEL_APP_GET_PANEL_API_BASE_PATH + panel_app_panel_id
    r = requests.get(url)
    json_data: Dict = r.json()
    # Panel App isn't very REST-ful - returns 200 for missing data but we'll return 404
    if detail := json_data.get("detail"):
        if detail == "Not found.":
            raise NotFound(detail=f"PanelApp couldn't find {panel_app_panel_id} ({url})")

    genes = json_data["genes"]
    name = json_data["name"]

    # Need to build evidence from panel app but using our gene names
    panel_app_gene_evidence = {}
    gene_names_list = []
    for gene_record in genes:
        gene_data = gene_record["gene_data"]
        gene_symbol = gene_data["gene_symbol"]
        panel_app_gene_evidence[gene_symbol] = gene_data
        gene_names_list.append(gene_symbol)

    gene_matcher = GeneSymbolMatcher()
    gene_list = FakeGeneList(name=name, user=None)
    gene_list_gene_symbols = gene_matcher.create_gene_list_gene_symbols(gene_list, gene_names_list, save=False)
    genelistgenesymbol_set = GeneListGeneSymbolSerializer(gene_list_gene_symbols, many=True).data  # sorted(glg_set)

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
            "absolute_url": url,
            "gene_evidence": gene_evidence}
    return data


def store_panel_app_panels_from_web(cached_web_resource):
    """ Used to build autocomplete list """

    # Now uses paging
    num_panels = 0
    url = PANEL_APP_BASE_URL + PANEL_APP_LIST_PANELS_PATH
    while url:
        logging.debug("Calling %s", url)
        r = requests.get(url)
        data = r.json()
        url = data["next"]

        for result in data["results"]:
            num_panels += 1
            pap = PanelAppPanel.objects.create(panel_id=result['id'],
                                               cached_web_resource=cached_web_resource,
                                               disease_group=result['disease_group'],
                                               disease_sub_group=result['disease_sub_group'],
                                               name=result['name'],
                                               current_version=result['version'])

            relevant_disorders = result['relevant_disorders']
            for rd in relevant_disorders:
                PanelAppPanelRelevantDisorders.objects.create(panel_app_panel=pap,
                                                              name=rd)

    cached_web_resource.description = f"{num_panels} panel app panels."
    cached_web_resource.save()
