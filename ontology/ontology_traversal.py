"""In-memory + DB-backed ontology graph traversal.

Both implementations share the BFS body in ``models_ontology.bfs_to_ontology``; only the
"give me adjacent relations" step differs. The DB version queries
``OntologyTermRelation``; the memory version walks adjacency dicts loaded
once.
"""
from collections import defaultdict
from typing import Protocol, Union

from genes.models import GeneSymbol
from ontology.models import (
    ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER,
    GeneDiseaseClassification,
    OntologyImportSource,
    OntologyRelation,
    OntologyRelationshipQualityFilter,
    OntologyService,
    OntologyTerm,
    OntologyTermRelation,
    OntologyVersion,
    PanelAppClassification,
)
from ontology.models.models_ontology import (
    EXCLUDED_SNAKE_RELATIONS,
    OntologySnake,
    OntologySnakes,
    bfs_to_ontology,
)


def _quality_filter_passes(otr: OntologyTermRelation,
                           gencc_ok: set[str], panelapp_ok: set[str]) -> bool:
    """Python mirror of OntologyRelationshipQualityFilter.filter_q for one row."""
    import_source = otr.from_import.import_source
    extra = otr.extra or {}
    strongest = extra.get("strongest_classification")
    if import_source == OntologyImportSource.GENCC:
        return strongest in gencc_ok
    if import_source == OntologyImportSource.PANEL_APP_AU:
        return strongest in panelapp_ok
    return True


class OntologyTraverser(Protocol):
    """Operations gene_annotation needs from an ontology graph for one
    ``OntologyVersion``. Two implementations: ``DbOntologyTraverser``
    (queries per call) and ``MemoryOntologyTraverser`` (adjacency dicts
    loaded once)."""

    def snake_from(self, term: OntologyTerm, to_ontology: OntologyService,
                   quality_filter: OntologyRelationshipQualityFilter = ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER,
                   max_depth: int = 1): ...

    def terms_for_gene_symbol(self, gene_symbol: Union[str, GeneSymbol],
                              desired_ontology: OntologyService, max_depth: int = 1,
                              quality_filter: OntologyRelationshipQualityFilter = ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER): ...

    def gene_disease_relations(self, gene_symbol: Union[str, GeneSymbol],
                               quality_filter: OntologyRelationshipQualityFilter = ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER) -> list[OntologyTermRelation]: ...


class DbOntologyTraverser:
    def __init__(self, ontology_version: OntologyVersion,
                 call_update_gene_relations: bool = True):
        self._ontology_version = ontology_version
        self._otr_qs = ontology_version.get_ontology_term_relations()
        self._call_update_gene_relations = call_update_gene_relations

    @property
    def ontology_version(self) -> OntologyVersion:
        return self._ontology_version

    def snake_from(self, term, to_ontology,
                   quality_filter=ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER,
                   max_depth=1):
        return OntologySnake.snake_from(term, to_ontology,
                                        quality_filter=quality_filter,
                                        max_depth=max_depth,
                                        otr_qs=self._otr_qs)

    def terms_for_gene_symbol(self, gene_symbol, desired_ontology, max_depth=1,
                              quality_filter=ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER):
        return self._ontology_version.terms_for_gene_symbol(
            gene_symbol, desired_ontology,
            max_depth=max_depth, quality_filter=quality_filter,
            call_update_gene_relations=self._call_update_gene_relations)

    def gene_disease_relations(self, gene_symbol,
                               quality_filter=ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER):
        return self._ontology_version.gene_disease_relations(
            gene_symbol, quality_filter=quality_filter,
            call_update_gene_relations=self._call_update_gene_relations)


class MemoryOntologyTraverser:
    """Loads one ``OntologyVersion``'s relation graph into adjacency dicts
    and serves the traverser API from it.

    Frozen snapshot — PanelApp live updates that the DB path triggers from
    ``OntologySnake.terms_for_gene_symbol`` are skipped here. Intended for
    batch runs that touch every HGNC term.
    """

    def __init__(self, ontology_version: OntologyVersion):
        self._ontology_version = ontology_version
        otr_qs = ontology_version.get_ontology_term_relations()
        edges = list(otr_qs.select_related("from_import", "source_term", "dest_term"))
        self._edges = edges
        self._by_source: dict[str, list[OntologyTermRelation]] = defaultdict(list)
        self._by_dest: dict[str, list[OntologyTermRelation]] = defaultdict(list)
        terms: dict[str, OntologyTerm] = {}
        for otr in edges:
            self._by_source[otr.source_term_id].append(otr)
            self._by_dest[otr.dest_term_id].append(otr)
            terms[otr.source_term_id] = otr.source_term
            terms[otr.dest_term_id] = otr.dest_term
        self._terms = terms
        self._hgnc_by_name = {
            t.name.upper(): t for t in terms.values()
            if t.ontology_service == OntologyService.HGNC and t.name
        }

    @property
    def ontology_version(self) -> OntologyVersion:
        return self._ontology_version

    @property
    def edge_count(self) -> int:
        return len(self._edges)

    def snake_from(self, term, to_ontology,
                   quality_filter=ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER,
                   max_depth=1):
        if term.ontology_service == to_ontology:
            return OntologySnakes([OntologySnake(source_term=term)])
        step_fn = self._make_step_fn(quality_filter)
        return bfs_to_ontology(term, to_ontology, max_depth, step_fn)

    def terms_for_gene_symbol(self, gene_symbol, desired_ontology, max_depth=1,
                              quality_filter=ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER):
        if isinstance(gene_symbol, GeneSymbol):
            symbol = gene_symbol.symbol
        else:
            symbol = str(gene_symbol)
        hgnc_term = self._hgnc_by_name.get(symbol.upper())
        if hgnc_term is None:
            return OntologySnakes([])
        return self.snake_from(hgnc_term, desired_ontology,
                               quality_filter=quality_filter, max_depth=max_depth)

    def gene_disease_relations(self, gene_symbol,
                               quality_filter=ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER):
        snakes = self.terms_for_gene_symbol(gene_symbol, OntologyService.MONDO,
                                            max_depth=0, quality_filter=quality_filter)
        return snakes.leaf_relations(ontology_relation=OntologyRelation.RELATED)

    def _make_step_fn(self, quality_filter: OntologyRelationshipQualityFilter):
        gencc_ok = GeneDiseaseClassification.get_above_min(quality_filter.min_gencc_strength)
        panelapp_ok = PanelAppClassification.get_above_min(quality_filter.min_panelapp_strength)

        def step(leaf_term_ids, seen_term_ids, exclude_hpo_endpoints):
            results: list[OntologyTermRelation] = []
            for leaf_id in leaf_term_ids:
                for otr in self._by_source.get(leaf_id, ()):
                    if otr.dest_term_id in seen_term_ids:
                        continue
                    if otr.relation in EXCLUDED_SNAKE_RELATIONS:
                        continue
                    if not _quality_filter_passes(otr, gencc_ok, panelapp_ok):
                        continue
                    if exclude_hpo_endpoints and otr.dest_term.ontology_service == OntologyService.HPO:
                        continue
                    results.append(otr)
                for otr in self._by_dest.get(leaf_id, ()):
                    if otr.source_term_id in seen_term_ids:
                        continue
                    if otr.relation in EXCLUDED_SNAKE_RELATIONS:
                        continue
                    if not _quality_filter_passes(otr, gencc_ok, panelapp_ok):
                        continue
                    if exclude_hpo_endpoints and otr.source_term.ontology_service == OntologyService.HPO:
                        continue
                    results.append(otr)
            return results

        return step


def get_ontology_traverser(ontology_version: OntologyVersion,
                           in_memory: bool = False,
                           call_update_gene_relations: bool = True) -> OntologyTraverser:
    if in_memory:
        # MemoryOntologyTraverser is a frozen snapshot — the live PanelApp hook
        # is unreachable from it regardless of the flag.
        return MemoryOntologyTraverser(ontology_version)
    return DbOntologyTraverser(ontology_version,
                               call_update_gene_relations=call_update_gene_relations)
