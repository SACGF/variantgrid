import re
from datetime import timedelta

from django.conf import settings
from typing import Union

from genes.models import GeneSymbol, PanelAppServer
from genes.panel_app import get_panel_app_results_by_gene_symbol_json, PANEL_APP_SEARCH_BY_GENES_BASE_PATH
from library.cache import timed_cache
from library.log_utils import report_exc_info
from library.utils import md5sum_str
from ontology.models import OntologyService, OntologyTerm, OntologyRelation, OntologyImportSource
from ontology.ontology_builder import OntologyBuilder, OntologyBuilderDataUpToDateException


def update_gene_relations(gene_symbol: Union[GeneSymbol, str]):
    if isinstance(gene_symbol, GeneSymbol):
        gene_symbol = gene_symbol.symbol
    return _update_gene_relations(gene_symbol)


@timed_cache(size_limit=2, ttl=10, quick_key_access=True)
def _update_gene_relations(gene_symbol: str):
    if not settings.PANEL_APP_CHECK_ENABLED:
        return

    # note that we only check PanelApp here, as other imports are done by file
    try:
        hgnc_term = OntologyTerm.get_gene_symbol(gene_symbol)
        panel_app = PanelAppServer.australia_instance()
        filename = panel_app.url + PANEL_APP_SEARCH_BY_GENES_BASE_PATH + gene_symbol
        ontology_builder = OntologyBuilder(filename=filename, context=str(gene_symbol), import_source=OntologyImportSource.PANEL_APP_AU, processor_version=2)
        try:
            ontology_builder.ensure_old(max_age=timedelta(days=settings.PANEL_APP_CACHE_DAYS))

            PANEL_APP_OMIM = re.compile(r"([0-9]{5,})")

            results = get_panel_app_results_by_gene_symbol_json(server=panel_app, gene_symbol=gene_symbol)
            response_hash = md5sum_str(str(results))

            ontology_builder.ensure_hash_changed(data_hash=response_hash)

            for panel_app_result in results:
                if evidence := panel_app_result.get('evidence'):
                    if "Expert Review Green" in evidence:  # only look at green panels
                        phenotype_row: str
                        for phenotype_row in panel_app_result.get("phenotypes", []):
                            # TODO look for terms other than OMIM in case panel app switches
                            for omim_match in PANEL_APP_OMIM.finditer(phenotype_row):
                                omim_int = int(omim_match[1])  # not always a valid id
                                omim_id = OntologyService.index_to_id(OntologyService.OMIM, omim_int)
                                if omim_term := OntologyTerm.objects.filter(id=omim_id).first():
                                    ontology_builder.add_ontology_relation(
                                        source_term_id=omim_term.id,
                                        dest_term_id=hgnc_term.id,
                                        relation=OntologyRelation.PANEL_APP_AU,
                                        extra={
                                            "phenotype_row": phenotype_row,
                                            "evidence": evidence
                                        }
                                    )

            ontology_builder.complete()

        except OntologyBuilderDataUpToDateException:
            pass
    except ValueError:
        report_exc_info()
        pass