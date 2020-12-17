import operator
from dataclasses import dataclass
from enum import Enum
from functools import reduce
from typing import Optional, List, Dict, Set

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


class GeneAssociationLevel(Enum):
    NONE = 0
    PARENT = 1
    DIRECT = 2

@dataclass
class ScorePart:
    category: str
    calculation: str
    score: float


class AutoMatchScore:

    def __init__(self, match: MonarchDiseaseOntology):
        self.match = match
        self.missing_word_count = 0.0
        self.search_word_count = 0.0
        self.superfluous_word_count = 0.0
        self.mondo_gene_association = GeneAssociationLevel.NONE

        self.score_parts = None
        self.score = None

    def apply_word_score(self, search_terms: Set[str]):
        match_terms = set(tokenize_condition_text(self.match.name))
        self.search_word_count = len(search_terms)

        self.superfluous_word_count = len(match_terms.difference(search_terms))
        self.missing_word_count = len(search_terms.difference(match_terms))

    def calculate(self):
        score = 10.0
        score_parts = list()
        missing_word_ratio = float(self.missing_word_count) / float(self.search_word_count)
        superfluous_word_ratio = float(self.superfluous_word_count) / float(self.search_word_count)

        if missing_word_ratio:
            score_parts.append(ScorePart(
                category="Missing word %",
                calculation=f"NEG Missing ratio : {missing_word_ratio} * 10",
                score=missing_word_ratio * -10.0
            ))

        if superfluous_word_ratio:
            score_parts.append(ScorePart(
                category="Superfluous word %",
                calculation=f"NEG Superfluous ratio : {superfluous_word_ratio} * 1",
                score=superfluous_word_ratio * -1.0
            ))

        if self.mondo_gene_association == GeneAssociationLevel.DIRECT:
            score_parts.append(ScorePart(
                category="Gene Association",
                calculation=f"Direction association",
                score=10.0
            ))
        if self.mondo_gene_association == GeneAssociationLevel.PARENT:
            score_parts.append(ScorePart(
                category="Gene Association",
                calculation=f"Child term has Gene association",
                score=5.0
            ))

        for score_part in score_parts:
            score += score_part.score

        self.score_parts = score_parts
        self.score = score

    def __lt__(self, other):
        return self.score < other.score

@dataclass
class MatchRequest:
    text: str
    gene_symbol: Optional[GeneSymbol]


def attempt_auto_match_for_text(request: MatchRequest) -> List[AutoMatchScore]:
    text = ConditionText.normalize(request.text)
    search_terms = set(tokenize_condition_text(text))

    if len(search_terms) >= 1:
        # this currently requires all words to be present
        qs = reduce(operator.and_, [Q(name__icontains=term) for term in search_terms])

        print(search_terms)

        matches_qs = MonarchDiseaseOntology.objects.filter(qs).order_by(Length('name').asc())  # find the smallest match that matches all the words
        match: MonarchDiseaseOntology

        scored_match_dict: Dict[MonarchDiseaseOntology, AutoMatchScore] = dict()

        for match in matches_qs[0:100]:
            scored_match_dict[match] = AutoMatchScore(match)

        if request.gene_symbol:
            for relationship in MonarchDiseaseOntologyGeneRelationship.objects.filter(gene_symbol=request.gene_symbol):
                gene_match = relationship.mondo
                scorer = AutoMatchScore(gene_match)
                scored_match_dict[gene_match] = scorer
                scorer.mondo_gene_association = GeneAssociationLevel.DIRECT

                for parent_relationship in MonarchDiseaseOntologyRelationship.objects.filter(subject=gene_match, relationship="is_a"):
                    gene_match_parent = parent_relationship.object
                    scorer = AutoMatchScore(gene_match)
                    scored_match_dict[gene_match] = scorer
                    scorer.mondo_gene_association = GeneAssociationLevel.PARENT

        values = scored_match_dict.values()
        for value in values:
            value.apply_word_score(search_terms)
            value.calculate()
        values = sorted(values, reverse=True)
        return values[0:3] # return up to top 3 results
    return []