import re
from datetime import timedelta
from typing import Union, Iterable

from genes.models import GeneSymbol, PanelAppServer
from genes.panel_app import get_panel_app_results_by_gene_symbol_json, PANEL_APP_SEARCH_BY_GENES_BASE_PATH
from library.log_utils import report_message
from library.utils import md5sum_str
from ontology.models import OntologySet, OntologyTermGeneRelation, OntologyTerm, OntologyTermRelation
from ontology.ontology_builder import OntologyBuilder, OntologyBuilderDataUpToDateException


def get_or_fetch_gene_relations(gene_symbol: Union[GeneSymbol, str]) -> Iterable[OntologyTermGeneRelation]:
    if not isinstance(gene_symbol, GeneSymbol):
        gene_symbol = GeneSymbol(symbol=gene_symbol)

    # note that we only check PanelApp here, as other imports are done by file

    panel_app = PanelAppServer.australia_instance()
    filename = panel_app.url + PANEL_APP_SEARCH_BY_GENES_BASE_PATH + gene_symbol.symbol

    ontology_builder = OntologyBuilder(filename=filename, context=str(gene_symbol), ontology_set=OntologySet.PANEL_APP_AU)
    try:
        ontology_builder.ensure_old(max_age=timedelta(days=30))

        PANEL_APP_OMIM = re.compile(r"([0-9]{5,})")

        results = get_panel_app_results_by_gene_symbol_json(server=panel_app, gene_symbol=gene_symbol)
        response_hash = md5sum_str(str(results))

        ontology_builder.ensure_hash_changed(data_hash=response_hash)

        for panel_app_result in results:
            phenotype_row: str
            for phenotype_row in panel_app_result.get("phenotypes", []):
                # TODO look for terms other than OMIM
                for omim_match in PANEL_APP_OMIM.finditer(phenotype_row):
                    found_mondo = False
                    omim_int = int(omim_match[1])  # not always a valid id
                    omim_id = OntologySet.index_to_id(OntologySet.OMIM, omim_int)
                    if omim_term := OntologyTerm.objects.filter(id=omim_id).first():
                        if mondo_term := OntologyTermRelation.mondo_version_of(omim_term):
                            # do we duplicate this for MONDO?
                            ontology_builder.add_gene_relation(
                                term_id=mondo_term.id,
                                gene_symbol_str=gene_symbol.symbol,
                                relation=OntologySet.PANEL_APP_AU,
                                extra={
                                    "phenotype_row": phenotype_row,
                                    "evidence": panel_app_result.get('evidence') or []
                                }
                            )
                            found_mondo = True
                    if not found_mondo:
                        report_message(f"Could not find a MONDO equiv of {omim_id} related to PanelApp {gene_symbol.symbol}")

        ontology_builder.complete()

    except OntologyBuilderDataUpToDateException:
        pass

    return OntologyTermGeneRelation.objects.filter(gene_symbol=gene_symbol)
