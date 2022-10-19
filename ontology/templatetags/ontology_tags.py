import itertools
import uuid
from dataclasses import dataclass
from typing import Optional, Union, Iterable, List, Iterator, Tuple
from uuid import UUID

from django.template import Library

from library.utils import segment
from ontology.models import OntologyTerm, OntologyTermRelation, GeneDiseaseClassification, OntologyService, \
    OntologySnake
from ontology.ontology_matching import OntologyMatch

register = Library()


@register.inclusion_tag("ontology/tags/ontology_match.html")
def ontology_meta(data: OntologyMatch):
    return {"ontology": data}


@register.inclusion_tag("ontology/tags/ontology_term.html")
def ontology_term(data: Union[OntologyTerm, str], show_link: bool = True, compact: bool = False):
    if isinstance(data, str):
        data = OntologyTerm.get_or_stub(data)

    is_gene = data.ontology_service == OntologyService.HGNC

    return {
        "term": data,
        "is_gene": is_gene,
        "show_link": show_link,
        "compact": compact
    }


@register.inclusion_tag("ontology/tags/ontology_relationship_table.html")
def ontology_relationship_table(relationships: OntologyTermRelation, reference_term: Optional[OntologyTerm] = None, other_term_title: str = "Other Term"):
    return {
        "other_term_title": other_term_title,
        "relationships": relationships,
        "reference_term": reference_term,
        "table_id": str(uuid.uuid4())
    }


@register.inclusion_tag("ontology/tags/ontology_relationship_row.html")
def ontology_relationship_row(relationship: OntologyTermRelation, reference_term: Optional[OntologyTerm] = None):
    low_quality = False
    quality: Optional[str] = None
    if extra := relationship.extra:
        if strongest := extra.get('strongest_classification'):
            allowed_set = GeneDiseaseClassification.get_above_min(GeneDiseaseClassification.STRONG)
            if strongest not in allowed_set:
                low_quality = True
                quality = strongest

    other_term = None
    if reference_term:
        other_term = relationship.other_end(reference_term)

    return {
        "other_term": other_term,
        "dest_term": other_term or relationship.dest_term,
        "relationship": relationship,
        "reversed": reference_term == relationship.dest_term,
        "low_quality": low_quality,
        "quality": quality,
        "reference_term": reference_term
    }


@dataclass
class GroupedSnakes:
    snakes: List[OntologySnake]
    destination: OntologyTerm


@dataclass
class GroupedSnakeRow:
    snake: OntologySnake
    row_span: int

    @staticmethod
    def yield_snakes(grouped_snakes: List['GroupedSnakes']) -> Iterator['GroupedSnakeRow']:
        for grouped_snake in grouped_snakes:
            is_first = True
            for snake in grouped_snake.snakes:
                row_span = len(grouped_snake.snakes) if is_first else 0
                is_first = False
                yield GroupedSnakeRow(snake=snake, row_span=row_span)


@register.inclusion_tag("ontology/tags/ontology_snake_table.html")
def ontology_snake_table(snakes: Iterable[OntologySnake], reference_term: Optional[OntologyTerm]):

    grouped: List[GroupedSnakes] = list()
    for leaf, snakes in itertools.groupby(snakes, lambda s: s.leaf_term):
        grouped.append(GroupedSnakes(snakes=list(snakes), destination=leaf))

    return {
        "table_id": str(uuid.uuid4()),
        "snakes": GroupedSnakeRow.yield_snakes(grouped),
        "reference_term": reference_term
    }


@register.inclusion_tag("ontology/tags/ontology_snake_row.html")
def ontology_snake_row(snake: OntologySnake, reference_term: Optional[OntologyTerm], row_span: int = 1):
    steps = snake.show_steps()
    source_term = snake.source_term
    dest_term = steps[-1].dest_term

    return {
        "source_term": source_term,
        "is_source_diff": (not reference_term) or (source_term != reference_term),
        "is_dest_diff": (not reference_term) or (dest_term != reference_term),
        "steps": steps,
        "dest_term": dest_term,
        "reference_term": reference_term,
        "row_span": row_span
    }
