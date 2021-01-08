import operator
import urllib
from dataclasses import dataclass
from enum import Enum
from functools import reduce
from typing import Optional, List, Dict, Set, Iterable

import requests
import re
from django.db.models import Q
from django.db.models.functions import Length
from lazy import lazy

from classification.models import ConditionText
from genes.models import GeneSymbol
from ontology.models import OntologyTermGeneRelation, OntologyTermRelation, OntologyTerm, OntologySet
from ontology.panel_app_ontology import get_or_fetch_gene_relations

SKIP_TERMS = {"", "disease", "disorder", "the", "a", "an", "and", "or", "for", "the", "type"}
SUB_TYPE = re.compile("[0-9]+[a-z]?")


def tokenize_condition_text(text: str):
    text = text.replace(";", " ").replace("-", " ").lower()
    tokens = [token.strip() for token in text.split(" ")]
    tokens = [token for token in tokens if token not in SKIP_TERMS]
    return tokens


class CheckAPI(Enum):
    NEVER = 0
    IF_NO_LOCAL = 1
    ALWAYS = 2


@dataclass
class MatchRequest:
    text: str
    gene_symbol: Optional[GeneSymbol]
    check_api: CheckAPI = CheckAPI.ALWAYS

    @lazy
    def tokenized_text(self) -> Set[str]:
        return set(tokenize_condition_text(ConditionText.normalize(self.text)))


class GeneAssociationLevel(Enum):
    MISSING = -1
    NONE = 0
    PARENT = 1
    DIRECT = 2


@dataclass
class ScorePart:
    calculation: str
    score: float


class AutoMatchScore:

    def __init__(self, match: OntologyTerm, request: MatchRequest):
        self.match = match
        self.request = request

        self.mondo_gene_association = GeneAssociationLevel.NONE
        self.mondo_gene_association_child: Optional[OntologyTerm] = None
        if request.gene_symbol:
            self.mondo_gene_association = GeneAssociationLevel.MISSING

        self.score_parts = None
        self.score = None

    def calculate(self):
        score = 0.0
        score_parts = list()

        match_terms = set(tokenize_condition_text(ConditionText.normalize(self.match.name))) - SKIP_TERMS
        superfluous_words = match_terms.difference(self.request.tokenized_text)
        missing_words = self.request.tokenized_text.difference(match_terms)

        missing_word_ratio = float(len(missing_words)) / float(len(match_terms))
        superfluous_word_ratio = float(len(superfluous_words)) / float(len(match_terms))

        def pretty_set(s: Set[str]) -> str:
            l = list(s)
            l = sorted(l)
            text = ", ".join(l)
            return f"\"{text}\""

        if self.mondo_gene_association == GeneAssociationLevel.NONE:
            score_parts.append(ScorePart(
                calculation=f"No gene symbol requested",
                score=100.0
            ))
        elif self.mondo_gene_association == GeneAssociationLevel.DIRECT:
            score_parts.append(ScorePart(
                calculation=f"Direction {self.request.gene_symbol.symbol} association",
                score=100.0
            ))
        elif self.mondo_gene_association == GeneAssociationLevel.PARENT:
            score_parts.append(ScorePart(
                calculation=f"Child term {self.mondo_gene_association_child.id} has {self.request.gene_symbol.symbol} association",
                score=75.0
            ))
        elif self.mondo_gene_association == GeneAssociationLevel.MISSING:
            score_parts.append(ScorePart(
                calculation=f"Missing association to {self.request.gene_symbol.symbol}",
                score=50.0
            ))

        if missing_word_ratio:
            score_parts.append(ScorePart(
                calculation=f"Missing words {pretty_set(missing_words)} @ ratio {missing_word_ratio:.2f} * 100",
                score=missing_word_ratio * -100.0
            ))

        if superfluous_word_ratio:
            score_parts.append(ScorePart(
                calculation=f"Superfluous words {pretty_set(superfluous_words)} @ ratio {superfluous_word_ratio:.2f} * 50",
                score=superfluous_word_ratio * -50.0
            ))

        if superfluous_words and not self.request.gene_symbol:
            for word in superfluous_words:
                if SUB_TYPE.match(word):
                    score_parts.append(ScorePart(
                        calculation=f"Superfluous words includes subtype \"{word}\"",
                        score=-20.0
                    ))

        for score_part in score_parts:
            score += score_part.score

        self.score_parts = score_parts
        self.score = score

    def __lt__(self, other):
        return self.score < other.score


class AutoMatchScores:

    def __init__(self, request: MatchRequest):
        self.request = request
        self.score_dict: Dict[OntologyTerm, AutoMatchScore] = dict()

    def __getitem__(self, item) -> AutoMatchScore:
        if existing := self.score_dict.get(item):
            return existing
        else:
            new_value = AutoMatchScore(item, self.request)
            self.score_dict[item] = new_value
            return new_value

    def values(self) -> Iterable[AutoMatchScore]:
        return self.score_dict.values()


def attempt_auto_match_for_text(request: MatchRequest) -> List[AutoMatchScore]:
    """
    TODO do we even want to entertain matching non MONDO for this?
    right now we only include it
    """
    search_terms = request.tokenized_text

    if not search_terms:
        return list()

    # this currently requires all words to be present
    qs = reduce(operator.and_, [Q(name__icontains=term) for term in search_terms])

    matches_qs = OntologyTerm.objects.filter(qs).filter(ontology_set=OntologySet.MONDO).order_by(Length('name').asc())  # find the smallest match that matches all the words
    match: OntologyTerm

    scored_match_dict = AutoMatchScores(request=request)

    had_local_matches = False
    for match in matches_qs[0:100]:
        scorer = scored_match_dict[match]  # just create it by referencing it, no data to add at this point
        had_local_matches = True

    if (not had_local_matches and request.check_api == CheckAPI.IF_NO_LOCAL) or request.check_api == CheckAPI.ALWAYS:
        url_text = ConditionText.normalize(request.text)
        url_text = urllib.parse.quote(url_text).replace('/', '%252F')
        results = requests.get(f'https://api.monarchinitiative.org/api/search/entity/autocomplete/{url_text}', {
            "prefix": "MONDO",
            "rows": 20,
            "minimal_tokenizer": "false",
            "category": "disease"
        }).json().get("docs")

        for result in results:
            o_id = result.get('id')
            label = result.get('label')
            if label:
                label = label[0]

            if api_match := OntologyTerm.objects.filter(id=o_id).first():
                scored_match_dict[api_match]

    if request.gene_symbol:
        relationship: OntologyTermGeneRelation
        for relationship in get_or_fetch_gene_relations(gene_symbol=request.gene_symbol):
            gene_match_term = relationship.term
            scorer = scored_match_dict[gene_match_term]
            scorer.mondo_gene_association = GeneAssociationLevel.DIRECT

            if gene_match_parent := OntologyTermRelation.parent_of(source_term=gene_match_term):
                scorer = scored_match_dict[gene_match_parent]
                scorer.mondo_gene_association = GeneAssociationLevel.PARENT
                scorer.mondo_gene_association_child = gene_match_term

    values = list(scored_match_dict.values())
    for value in values:
        value.calculate()
    values = sorted(values, reverse=True)
    return values[0:10] # return up to top 3 results