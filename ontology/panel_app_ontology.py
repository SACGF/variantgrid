import re
from datetime import timedelta
from typing import Union, Iterable, Optional

from genes.models import GeneSymbol, PanelAppServer, HGNCGeneNames
from genes.panel_app import get_panel_app_results_by_gene_symbol_json, PANEL_APP_SEARCH_BY_GENES_BASE_PATH
from library.utils import md5sum_str
from ontology.models import OntologyService, OntologyTerm, OntologyRelation
from ontology.ontology_builder import OntologyBuilder, OntologyBuilderDataUpToDateException


def update_gene_relations(gene_symbol: Union[GeneSymbol, str]):
    if isinstance(gene_symbol, GeneSymbol):
        gene_symbol = gene_symbol.symbol

    # note that we only check PanelApp here, as other imports are done by file

    panel_app = PanelAppServer.australia_instance()
    filename = panel_app.url + PANEL_APP_SEARCH_BY_GENES_BASE_PATH + gene_symbol
    hgnc_term = OntologyTerm.get_gene_symbol(gene_symbol)

    ontology_builder = OntologyBuilder(filename=filename, context=str(gene_symbol), ontology_service=OntologyService.PANEL_APP_AU)
    try:
        ontology_builder.ensure_old(max_age=timedelta(days=30))

        PANEL_APP_OMIM = re.compile(r"([0-9]{5,})")

        results = get_panel_app_results_by_gene_symbol_json(server=panel_app, gene_symbol=gene_symbol)
        response_hash = md5sum_str(str(results))

        ontology_builder.ensure_hash_changed(data_hash=response_hash)

        for panel_app_result in results:
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
                                "evidence": panel_app_result.get('evidence') or []
                            }
                        )

        ontology_builder.complete()

    except OntologyBuilderDataUpToDateException:
        pass
