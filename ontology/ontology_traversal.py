"""In-memory + DB-backed ontology graph traversal.

Both implementations share the BFS body in ``bfs_to_ontology``; only the
"give me adjacent relations" step differs. The DB version queries
``OntologyTermRelation``; the memory version walks adjacency dicts loaded
once.
"""
from collections import defaultdict
from typing import Callable, Optional, Protocol, Union

from django.db.models import Q

from genes.models import GeneSymbol
from ontology.models import (
    GeneDiseaseClassification,
    OntologyImportSource,
    OntologyRelation,
    OntologyRelationshipQualityFilter,
    OntologyService,
    OntologyTerm,
    OntologyTermRelation,
    OntologyVersion,
    PanelAppClassification,
    ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER,
)


_EXCLUDED_RELATIONS = frozenset({
    OntologyRelation.IS_A,
    OntologyRelation.EXACT_SYNONYM,
    OntologyRelation.RELATED_SYNONYM,
    OntologyRelation.XREF,
})


def bfs_to_ontology(
    start_term: OntologyTerm,
    to_ontology: OntologyService,
    max_depth: int,
    step_fn: Callable[[set[str], set[str], bool], list[OntologyTermRelation]],
):
    """BFS from ``start_term`` toward ``to_ontology``.

    ``step_fn(leaf_term_ids, seen_term_ids, exclude_hpo_endpoints)`` returns
    the adjacent ``OntologyTermRelation`` rows to expand. Quality filtering,
    the excluded-relation set, and HPO-endpoint exclusion are baked into
    ``step_fn``.
    """
    from ontology.models import OntologySnake, OntologySnakes  # avoid circular import

    if start_term.ontology_service == to_ontology:
        return OntologySnakes([OntologySnake(source_term=start_term)])

    seen: set[str] = {start_term.id}
    new_snakes: list = [OntologySnake(source_term=start_term)]
    valid_snakes: list = []
    exclude_hpo_endpoints = (to_ontology == OntologyService.HGNC)

    while new_snakes:
        snakes = list(new_snakes)
        new_snakes = []
        by_leafs: dict[str, OntologySnake] = {}
        for snake in snakes:
            existing = by_leafs.get(snake.leaf_term.id)
            if existing is None or len(snake.paths) < len(existing.paths):
                by_leafs[snake.leaf_term.id] = snake

        all_relations = step_fn(set(by_leafs.keys()), set(seen), exclude_hpo_endpoints)

        for relation in all_relations:
            snake = by_leafs.get(relation.source_term_id) or by_leafs.get(relation.dest_term_id)
            if snake is None:
                continue
            if snake.leaf_term.id not in (relation.source_term_id, relation.dest_term_id):
                continue
            other_term = relation.other_end(snake.leaf_term)

            ontology_services = {snake.leaf_term.ontology_service, other_term.ontology_service}
            if (OntologyService.MONDO in ontology_services
                    and OntologyService.OMIM in ontology_services
                    and relation.relation != OntologyRelation.EXACT):
                continue

            new_snake = snake.snake_step(relation)
            if other_term.ontology_service == to_ontology:
                valid_snakes.append(new_snake)
                continue
            if len(new_snake.paths) <= max_depth:
                new_snakes.append(new_snake)
            seen.add(other_term.id)

    valid_snakes.sort()
    return OntologySnakes(valid_snakes)


def _make_db_step_fn(otr_qs, quality_filter: OntologyRelationshipQualityFilter):
    q_relation = ~Q(relation__in=_EXCLUDED_RELATIONS) & quality_filter.filter_q

    def step(leaf_term_ids, seen_term_ids, exclude_hpo_endpoints):
        outgoing = otr_qs.filter(source_term_id__in=leaf_term_ids).exclude(dest_term_id__in=seen_term_ids).filter(q_relation)
        incoming = otr_qs.filter(dest_term_id__in=leaf_term_ids).exclude(source_term_id__in=seen_term_ids).filter(q_relation)
        if exclude_hpo_endpoints:
            outgoing = outgoing.exclude(dest_term__ontology_service=OntologyService.HPO)
            incoming = incoming.exclude(source_term__ontology_service=OntologyService.HPO)
        return list(outgoing) + list(incoming)

    return step


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
        from ontology.models import OntologySnake
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
        from ontology.models import OntologySnake, OntologySnakes
        if term.ontology_service == to_ontology:
            return OntologySnakes([OntologySnake(source_term=term)])
        step_fn = self._make_step_fn(quality_filter)
        return bfs_to_ontology(term, to_ontology, max_depth, step_fn)

    def terms_for_gene_symbol(self, gene_symbol, desired_ontology, max_depth=1,
                              quality_filter=ONTOLOGY_RELATIONSHIP_STANDARD_QUALITY_FILTER):
        from ontology.models import OntologySnakes
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
                    if otr.relation in _EXCLUDED_RELATIONS:
                        continue
                    if not _quality_filter_passes(otr, gencc_ok, panelapp_ok):
                        continue
                    if exclude_hpo_endpoints and otr.dest_term.ontology_service == OntologyService.HPO:
                        continue
                    results.append(otr)
                for otr in self._by_dest.get(leaf_id, ()):
                    if otr.source_term_id in seen_term_ids:
                        continue
                    if otr.relation in _EXCLUDED_RELATIONS:
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
