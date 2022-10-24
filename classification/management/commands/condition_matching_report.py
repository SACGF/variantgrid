from django.core.management import BaseCommand

from library.utils import ExportRow, export_column, delimited_row
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional
from classification.models import ConditionTextMatch
from genes.models import GeneSymbol
from ontology.models import OntologyTerm, OntologySnake, OntologyService


@dataclass
class MatchKey:
    ontology_term: OntologyTerm
    gene_symbol: GeneSymbol

    def __hash__(self):
        return hash(self.ontology_term.id + "/" + self.gene_symbol.symbol)


class MatchLevel(int, Enum):
    NO_MATCH = 0
    NON_EXACT = 1
    EXACT = 2


@dataclass
class MatchRelationship:
    level: MatchLevel
    snake: OntologySnake

    def __lt__(self, other):
        return self.level < other.level


@dataclass
class MatchValue(ExportRow):
    key: MatchKey
    matches: List[ConditionTextMatch] = field(default_factory=list)
    match_level: Optional[MatchLevel] = MatchLevel.NO_MATCH
    matching_relationship: Optional[OntologySnake] = None

    def __lt__(self, other):
        return self.key.ontology_term < other.key.ontology_term

    @export_column("Plain Text")
    def export_texts(self):
        return ";".join(f"{match.condition_text.normalized_text} " for match in self.matches)

    @export_column("Labs")
    def export_labs(self):
        return ";".join(f"{match.condition_text.lab.name} " for match in self.matches)

    @export_column("Users")
    def export_users(self):
        return ";".join(f"{match.last_edited_by} " for match in self.matches)

    @export_column("Term")
    def export_term(self):
        return self.key.ontology_term.id

    @export_column("Name")
    def export_name(self):
        return self.key.ontology_term.name

    @export_column("Gene Symbol")
    def export_gene_symbol(self):
        return self.key.gene_symbol.symbol

    @export_column("Match Level")
    def export_match_level(self):
        return MatchLevel(self.match_level).name

    def __str__(self):
        return f"{self.key.ontology_term.id}\t{self.key.ontology_term.name}\t{self.key.gene_symbol.symbol}\t{MatchLevel(self.match_level).name}\t{self.matching_relationship}"


def is_bad_snake(snake: OntologySnake):
    for step in snake.show_steps():
        if step.source_term.ontology_service == OntologyService.MONDO and step.dest_term.ontology_service == OntologyService.OMIM and step.relation.relation != "exact":
            return True
    return False


def get_last_term_symbol_relationships(ontology_term: OntologyTerm) -> Dict[GeneSymbol, MatchRelationship]:
    data: Dict[GeneSymbol, MatchRelationship] = dict()
    for snake in OntologySnake.snake_from(ontology_term, to_ontology=OntologyService.HGNC):
        gene_symbol = GeneSymbol(snake.leaf_term.name)
        match_level = MatchLevel.EXACT
        if is_bad_snake(snake):
            match_level = MatchLevel.NON_EXACT

        existing_match_level = MatchLevel.NO_MATCH
        if existing_match_relation := data.get(gene_symbol):
            existing_match_level = existing_match_relation.level

        if match_level >= existing_match_level:
            data[gene_symbol] = MatchRelationship(level=match_level, snake=snake)

    return data


class Command(BaseCommand):

    def handle(self, *args, **options):

        all_matches: Dict[MatchKey, MatchValue] = dict()

        for ct in ConditionTextMatch.objects.filter(gene_symbol__isnull=False, mode_of_inheritance__isnull=True, classification__isnull=True):
            if conditions := ct.condition_xref_terms:
                if len(conditions) == 1:
                    match_key = MatchKey(ontology_term=conditions[0], gene_symbol=GeneSymbol.objects.get(symbol=ct.gene_symbol))
                    match_value: MatchValue
                    if existing := all_matches.get(match_key):
                        match_value = existing
                    else:
                        match_value = MatchValue(key=match_key)
                    match_value.matches.append(ct)
                    all_matches[match_key] = match_value

        all_values = sorted(all_matches.values())

        last_term: Optional[OntologyTerm] = None
        last_term_symbol_relationships: Dict[GeneSymbol, MatchRelationship] = dict()

        print(delimited_row(MatchValue.csv_header(), delimiter="\t"), end='')
        for match_value in all_values:
            this_term = match_value.key.ontology_term
            if this_term != last_term:
                last_term_symbol_relationships = get_last_term_symbol_relationships(this_term)
                last_term = this_term

            if match_relationship := last_term_symbol_relationships.get(match_value.key.gene_symbol):
                match_value.match_level = match_relationship.level
                match_value.matching_relationship = match_relationship.snake

            if match_value.match_level != MatchLevel.EXACT:
                print(delimited_row(match_value.to_csv(), delimiter="\t"), end='')
