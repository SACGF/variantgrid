import logging
from typing import Dict, Optional

import requests
from rest_framework.exceptions import NotFound

from annotation.models import CachedWebResource
from genes.gene_matching import GeneSymbolMatcher
from genes.models import PanelAppPanelRelevantDisorders, PanelAppPanel, PanelAppServer, PanelAppPanelLocalCache, \
    PanelAppPanelLocalCacheGeneSymbol, GeneSymbol, create_fake_gene_list
from genes.serializers import GeneListGeneSymbolSerializer
from library.constants import MINUTE_SECS

PANEL_APP_PREFIX = "panel-app"
PANEL_APP_LIST_PANELS_PATH = "/api/v1/panels/"
PANEL_APP_SEARCH_BY_GENES_BASE_PATH = "/api/v1/genes/"


def get_panel_app_results_by_gene_symbol_json(server: PanelAppServer, gene_symbol) -> Optional[Dict]:
    url = server.url + PANEL_APP_SEARCH_BY_GENES_BASE_PATH + str(gene_symbol)
    r = requests.get(url, timeout=MINUTE_SECS)
    results = None
    if r.ok:
        data = r.json()
        results = data.get("results")
    r.close()
    return results


def _get_panel_app_panel_api_json(panel_app_panel):
    r = requests.get(panel_app_panel.url, timeout=MINUTE_SECS)
    json_data: Dict = r.json()
    # Panel App isn't very REST-ful - returns 200 for missing data but we'll return 404
    if detail := json_data.get("detail"):
        if detail == "Not found.":
            raise NotFound(detail=f"PanelApp couldn't find {panel_app_panel.panel_id} ({url})")

    return json_data


def get_panel_app_panel_as_gene_list_json(panel_app_panel_id):
    panel_app_panel = PanelAppPanel.objects.get(pk=panel_app_panel_id)
    json_data = _get_panel_app_panel_api_json(panel_app_panel)

    genes = json_data["genes"]
    name = json_data["name"]

    # Need to build evidence from panel app but using our gene names
    panel_app_gene_evidence = {}
    gene_names_list = []
    for gene_record in genes:
        gene_symbol = gene_record["gene_data"]["gene_symbol"]
        panel_app_gene_evidence[gene_symbol] = gene_record
        gene_names_list.append(gene_symbol)

    gene_matcher = GeneSymbolMatcher()
    gene_list = create_fake_gene_list(name=name, user=None)
    gene_list_gene_symbols = gene_matcher.create_gene_list_gene_symbols(gene_list, gene_names_list, save=False)
    genelistgenesymbol_set = GeneListGeneSymbolSerializer(gene_list_gene_symbols, many=True).data  # sorted(glg_set)

    gene_evidence = {}
    for glgs in gene_list_gene_symbols:
        if glgs.gene_symbol_id:
            evidence = panel_app_gene_evidence[glgs.original_name]
            gene_evidence[glgs.gene_symbol_id] = evidence
            # TODO: Handle unmatched symbols??

    data = {
        "pk": f"{PANEL_APP_PREFIX}-{panel_app_panel_id}",
        "category": {"name": "PanelApp", "icon_css_class": panel_app_panel.server.icon_css_class},
        "name": name,
        "import_status": "S",
        "genelistgenesymbol_set": genelistgenesymbol_set,
        "can_write": False,
        "absolute_url": panel_app_panel.url,
        "gene_evidence": gene_evidence,
    }
    return data


def _get_or_update_panel_app_panel(server, json_data):
    pap, _ = PanelAppPanel.objects.update_or_create(server=server,
                                                    panel_id=json_data['id'],
                                                    defaults={
                                                        "disease_group": json_data['disease_group'],
                                                        "disease_sub_group": json_data['disease_sub_group'],
                                                        "name": json_data['name'],
                                                        "status": json_data['status'],
                                                        "current_version": json_data['version'],
                                                    })

    relevant_disorders = json_data['relevant_disorders']
    pap.panelapppanelrelevantdisorders_set.exclude(name__in=relevant_disorders).delete()
    for rd in relevant_disorders:
        PanelAppPanelRelevantDisorders.objects.get_or_create(panel_app_panel=pap, name=rd)

    return pap


def store_panel_app_panels_from_web(server: PanelAppServer, cached_web_resource: CachedWebResource):
    """ Used to build autocomplete list """

    # Now uses paging
    num_panels = 0
    url = server.url + PANEL_APP_LIST_PANELS_PATH
    while url:
        logging.debug("Calling %s", url)
        r = requests.get(url, timeout=MINUTE_SECS)
        data = r.json()
        url = data["next"]

        for result in data["results"]:
            num_panels += 1
            _get_or_update_panel_app_panel(server, result)

    cached_web_resource.description = f"{num_panels} panel app panels."
    cached_web_resource.save()


def get_panel_app_local_cache(panel_app_panel: PanelAppPanel) -> PanelAppPanelLocalCache:
    """ Gets or creates local cache of a panel app panel so it can be used as a GeneList """

    # Attempt to use cache if recent and present, otherwise fall through and do a query
    try:
        existing = PanelAppPanelLocalCache.objects.get(panel_app_panel=panel_app_panel,
                                                       version=panel_app_panel.current_version)
        if panel_app_panel.cache_valid:
            return existing
    except PanelAppPanelLocalCache.DoesNotExist:
        existing = None

    json_data = _get_panel_app_panel_api_json(panel_app_panel)
    genes = json_data["genes"]
    version = json_data["version"]

    # Update panel_app_panel to latest (even if just bumping modified date)
    _get_or_update_panel_app_panel(panel_app_panel.server, json_data)

    if existing and existing.version == version:
        return existing  # cache still valid

    pap_lc = PanelAppPanelLocalCache.objects.create(panel_app_panel=panel_app_panel,
                                                    version=version)
    existing_uc_symbols = GeneSymbol.get_upper_case_lookup()
    new_symbols = []
    pap_lc_genes = []
    for api_record in genes:
        gene_symbol = api_record["gene_data"]["gene_symbol"]
        if gene_symbol.upper() not in existing_uc_symbols:
            new_symbols.append(gene_symbol)
        record = PanelAppPanelLocalCacheGeneSymbol(panel_app_local_cache=pap_lc,
                                                   gene_symbol_id=gene_symbol,
                                                   data=api_record)
        pap_lc_genes.append(record)

    if new_symbols:
        GeneSymbol.objects.bulk_create(new_symbols, ignore_conflicts=True, batch_size=2000)
    if pap_lc_genes:
        PanelAppPanelLocalCacheGeneSymbol.objects.bulk_create(pap_lc_genes, batch_size=2000)
    return pap_lc
