import operator
from dataclasses import dataclass
from enum import Enum
from functools import reduce
from typing import Optional, List, Dict, Set, Iterable

from django.db.models import Q
from django.db.models.functions import Length
from lazy import lazy

from annotation.models import MonarchDiseaseOntology, MonarchDiseaseOntologyGeneRelationship, \
    MonarchDiseaseOntologyRelationship
from classification.models import ConditionTextMatch, ConditionText
from genes.models import GeneSymbol


SKIP_TERMS = {"", "disease", "disorder", "the", "a", "an", "for", "the"}

def tokenize_condition_text(text: str):
    text = text.replace(";", " ").replace("-", " ").lower()
    tokens = [token.strip() for token in text.split(" ")]
    tokens = [token for token in tokens if token not in SKIP_TERMS]
    return tokens


@dataclass
class MatchRequest:
    text: str
    gene_symbol: Optional[GeneSymbol]

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

    def __init__(self, match: MonarchDiseaseOntology, request: MatchRequest):
        self.match = match
        self.request = request

        self.mondo_gene_association = GeneAssociationLevel.NONE
        self.mondo_gene_association_child: Optional[MonarchDiseaseOntology] = None
        if request.gene_symbol:
            self.mondo_gene_association = GeneAssociationLevel.MISSING

        self.score_parts = None
        self.score = None


    def calculate(self):
        score = 0.0
        score_parts = list()

        match_terms = set(tokenize_condition_text(ConditionText.normalize(self.match.name)))
        superfluous_words = match_terms.difference(self.request.tokenized_text)
        missing_words = self.request.tokenized_text.difference(match_terms)

        missing_word_ratio = float(len(missing_words)) / float(len(match_terms))
        superfluous_word_ratio = float(len(superfluous_words)) / float(len(match_terms))

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
                calculation=f"Child term {self.mondo_gene_association_child.id_str} has {self.request.gene_symbol.symbol} association",
                score=75.0
            ))
        elif self.mondo_gene_association == GeneAssociationLevel.MISSING:
            score_parts.append(ScorePart(
                calculation=f"Missing association to {self.request.gene_symbol.symbol}",
                score=50.0
            ))

        if missing_word_ratio:
            score_parts.append(ScorePart(
                calculation=f"Missing words {missing_words} @ ratio {missing_word_ratio:.2f} * 100",
                score=missing_word_ratio * -100.0
            ))

        if superfluous_word_ratio:
            score_parts.append(ScorePart(
                calculation=f"Superfluous words {superfluous_words} @ ratio {superfluous_word_ratio:.2f} * 10",
                score=superfluous_word_ratio * -10.0
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
        self.score_dict: Dict[MonarchDiseaseOntology, AutoMatchScore] = dict()

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
    search_terms = request.tokenized_text

    if not search_terms:
        return list()

    # this currently requires all words to be present
    qs = reduce(operator.and_, [Q(name__icontains=term) for term in search_terms])

    matches_qs = MonarchDiseaseOntology.objects.filter(qs).order_by(Length('name').asc())  # find the smallest match that matches all the words
    match: MonarchDiseaseOntology

    scored_match_dict = AutoMatchScores(request=request)

    for match in matches_qs[0:100]:
        scorer = scored_match_dict[match]  # just create it by referencing it, no data to add at this point

    if request.gene_symbol:
        for relationship in MonarchDiseaseOntologyGeneRelationship.objects.filter(gene_symbol=request.gene_symbol):
            gene_match = relationship.mondo
            scorer = scored_match_dict[gene_match]
            scorer.mondo_gene_association = GeneAssociationLevel.DIRECT

            for parent_relationship in MonarchDiseaseOntologyRelationship.objects.filter(subject=gene_match, relationship="is_a"):
                gene_match_parent = parent_relationship.object
                scorer = scored_match_dict[gene_match_parent]
                scorer.mondo_gene_association = GeneAssociationLevel.PARENT
                scorer.mondo_gene_association_child = gene_match

    values = list(scored_match_dict.values())
    for value in values:
        value.calculate()
    values = sorted(values, reverse=True)
    return values[0:5] # return up to top 3 results