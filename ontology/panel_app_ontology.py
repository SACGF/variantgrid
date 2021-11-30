from datetime import timedelta
from typing import Union

from django.conf import settings

from genes.models import GeneSymbol, PanelAppServer
from genes.panel_app import get_panel_app_results_by_gene_symbol_json, PANEL_APP_SEARCH_BY_GENES_BASE_PATH
from library.cache import timed_cache
from library.log_utils import report_exc_info, report_message
from library.utils import md5sum_str
from ontology.models import OntologyTerm, OntologyRelation, OntologyImportSource
from ontology.ontology_builder import OntologyBuilder, OntologyBuilderDataUpToDateException

# increment if you change the logic of parsing ontology terms from PanelApp
# which will then effectively nullify the cache so the new logic is run
PANEL_APP_API_PROCESSOR_VERSION = 4


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
        ontology_builder = OntologyBuilder(filename=filename, context=str(gene_symbol),
                                           import_source=OntologyImportSource.PANEL_APP_AU,
                                           processor_version=PANEL_APP_API_PROCESSOR_VERSION)
        try:
            ontology_builder.ensure_old(max_age=timedelta(days=settings.PANEL_APP_CACHE_DAYS))

            results = get_panel_app_results_by_gene_symbol_json(server=panel_app, gene_symbol=gene_symbol)
            response_hash = md5sum_str(str(results))

            ontology_builder.ensure_hash_changed(data_hash=response_hash)

            for panel_app_result in results:
                if evidence := panel_app_result.get('evidence'):
                    if "Expert Review Green" in evidence:  # only look at green panels
                        phenotype_row: str
                        for phenotype_row in panel_app_result.get("phenotypes", []):

                            def add_term_if_valid(full_id: str):
                                nonlocal ontology_builder
                                nonlocal hgnc_term
                                if term := OntologyTerm.objects.filter(id=full_id).first():
                                    ontology_builder.add_ontology_relation(
                                        source_term_id=term.id,
                                        dest_term_id=hgnc_term.id,
                                        relation=OntologyRelation.PANEL_APP_AU,
                                        extra={
                                            "phenotype_row": phenotype_row,
                                            "evidence": evidence
                                        })
                                else:
                                    report_message("Found ontology term from PanelApp not in DB", level="error", extra_data={"target": full_id, "gene_symbol": str(gene_symbol)})

                            from annotation.regexes import db_ref_regexes, DbRegexes
                            for result in db_ref_regexes.search(phenotype_row, default_regex=DbRegexes.OMIM):
                                if result.cregx in (DbRegexes.OMIM, DbRegexes.MONDO):
                                    add_term_if_valid(result.id_fixed)

            ontology_builder.complete()

        except OntologyBuilderDataUpToDateException:
            pass
    except ValueError:
        report_exc_info()
