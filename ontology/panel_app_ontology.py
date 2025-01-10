import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import timedelta
from typing import Union, Optional, Any

from django.conf import settings
from django.db import transaction

from genes.models import GeneSymbol, PanelAppServer
from genes.panel_app import get_panel_app_results_by_gene_symbol_json, PANEL_APP_SEARCH_BY_GENES_BASE_PATH
from library.cache import timed_cache
from library.log_utils import report_exc_info, report_message
from library.utils import md5sum_str
from ontology.models import OntologyTerm, OntologyRelation, OntologyImportSource, OntologyTermRelation, \
    OntologyTermStatus, OntologyIdNormalized, PanelAppClassification
from ontology.ontology_builder import OntologyBuilder, OntologyBuilderDataUpToDateException

# increment if you change the logic of parsing ontology terms from PanelApp
# which will then effectively nullify the cache so the new logic is run
PANEL_APP_API_PROCESSOR_VERSION = 11
# with look ahead and behind to make sure we're not in a 7-digit number
ABANDONED_OMIM_RE = re.compile('(?<![0-9])([0-9]{6})(?![0-9])')


def update_gene_relations(gene_symbol: Union[GeneSymbol, str]):
    if isinstance(gene_symbol, GeneSymbol):
        gene_symbol = gene_symbol.symbol
    return _update_gene_relations(gene_symbol)


@dataclass(frozen=True)
class PanelAppResult:

    ontology_ids: set[str]
    max_panel_app_strength: Optional[PanelAppClassification]
    raw_data: dict
    hash_str: str

    @staticmethod
    def parse_data(raw_data: dict) -> 'PanelAppResult':
        phenotypes = raw_data.get('phenotypes', [])
        evidences = raw_data.get('evidence', [])
        hash_str = ""

        all_terms = set()
        for phenotype_row in phenotypes:
            hash_str += phenotype_row + ";"
            found_term = False
            from annotation.regexes import db_ref_regexes, DbRegexes
            for result in db_ref_regexes.search(phenotype_row):
                if result.cregx in (DbRegexes.OMIM, DbRegexes.MONDO):
                    all_terms.add(result.id_fixed)
            if not found_term:
                # just look for abandoned 6-digit numbers
                for omim_id in ABANDONED_OMIM_RE.finditer(phenotype_row):
                    all_terms.add(f"OMIM:{omim_id.group(1)}")

        max_strength: Optional[PanelAppClassification] = None
        for evidence in evidences:
            hash_str += evidence + ";"
            try:
                panel_app_strength = PanelAppClassification.get_by_label_pac(evidence)
                if max_strength is None or max_strength < panel_app_strength:
                    max_strength = panel_app_strength
            except ValueError:
                # not actually a PanelAppStrength
                pass

        return PanelAppResult(
            ontology_ids=all_terms,
            max_panel_app_strength=max_strength,
            raw_data={
                # this is just a subset of the full data
                "phenotypes": phenotypes,
                "evidence": evidences
            },
            hash_str=hash_str
        )


@timed_cache(size_limit=2, ttl=10, quick_key_access=True)
def _update_gene_relations(gene_symbol: str):
    if not settings.GENE_RELATION_PANEL_APP_LIVE_UPDATE:
        return

    def is_empty_results(data: Any):
        return isinstance(data, list) and len(data) == 0

    # note that we only check PanelApp here, as other imports are done by file
    try:
        hgnc_term = OntologyTerm.get_gene_symbol(gene_symbol)
        panel_app = PanelAppServer.australia_instance()
        filename = panel_app.url + PANEL_APP_SEARCH_BY_GENES_BASE_PATH + gene_symbol
        ontology_builder = OntologyBuilder(filename=filename, context=str(gene_symbol),
                                           import_source=OntologyImportSource.PANEL_APP_AU,
                                           processor_version=PANEL_APP_API_PROCESSOR_VERSION,
                                           versioned=False)
        try:
            ontology_builder.ensure_old(max_age=timedelta(days=settings.PANEL_APP_CACHE_DAYS))
            if ontology_builder.versioned:
                raise ValueError("Can't do PanelAppAU with a versioned OntologyBuilder")

            alias_symbol: Optional[str] = None
            results_json = get_panel_app_results_by_gene_symbol_json(server=panel_app, gene_symbol=gene_symbol)
            if is_empty_results(results_json):
                # try aliases see if any of those work
                try:
                    gene_symbol_obj = GeneSymbol.objects.get(symbol=gene_symbol)
                    for alias in gene_symbol_obj.alias_meta.aliases_in:
                        if alias.other_symbol_in_database and not alias.different_genes:
                            alias_results_json = get_panel_app_results_by_gene_symbol_json(server=panel_app,
                                                                                     gene_symbol=alias.other_symbol)
                            if not is_empty_results(alias_results_json):
                                alias_symbol = alias.other_symbol
                                report_message(message=f"PanelAppAU no results for {gene_symbol}, substituting alias {alias.other_symbol}")
                                results_json = alias_results_json
                                break
                except GeneSymbol.DoesNotExist:
                    pass

            response_hash = md5sum_str(str(results_json))
            ontology_builder.ensure_hash_changed(data_hash=response_hash)

            OntologyTermRelation.objects.filter(
                dest_term_id=hgnc_term.id,
                relation=OntologyRelation.PANEL_APP_AU
            ).delete()

            by_ontology_id = defaultdict(list)

            for panel_app_result_json in results_json:
                parsed_result = PanelAppResult.parse_data(panel_app_result_json)
                if parsed_result.max_panel_app_strength and parsed_result.ontology_ids:
                    for ontology_id in parsed_result.ontology_ids:
                        by_ontology_id[ontology_id].append(parsed_result)

            with transaction.atomic():
                for ontology_id, parsed_results in by_ontology_id.items():
                    full_id = OntologyIdNormalized.normalize(ontology_id).full_id

                    term, created = ontology_builder.add_term(
                        term_id=full_id,
                        name="Unknown Term",
                        primary_source=False,
                        trusted_source=False
                    )
                    if created and term.status == OntologyTermStatus.STUB:
                        report_message("Found ontology term from PanelApp not in DB", level="error",
                                       extra_data={"target": full_id, "gene_symbol": str(gene_symbol)})

                    max_strength: Optional[PanelAppClassification] = None
                    all_data = []
                    unique_raw_data = set()
                    for parsed_result in parsed_results:
                        if max_strength is None or max_strength < parsed_result.max_panel_app_strength:
                            max_strength = parsed_result.max_panel_app_strength
                        if parsed_result.hash_str not in unique_raw_data:
                            all_data.append(parsed_result.raw_data)
                            unique_raw_data.add(parsed_result.hash_str)

                    extra = {
                        "strongest_classification": max_strength.label,
                        "phenotypes_and_evidence": all_data
                    }
                    if alias_symbol:
                        extra["using_alias"] = alias_symbol

                    ontology_builder.add_ontology_relation(
                        source_term_id=term.id,
                        dest_term_id=hgnc_term.id,
                        relation=OntologyRelation.PANEL_APP_AU,
                        extra=extra
                    )
                ontology_builder.complete(verbose=False)

        except OntologyBuilderDataUpToDateException:
            pass
    except ValueError:
        report_exc_info()
