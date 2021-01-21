import operator
import re
from dataclasses import dataclass
from functools import reduce
from typing import Dict, Optional, List, Any, Set

import requests
from django.db.models import Q
from django.db.models.functions import Length
from django.urls import reverse
from lazy import lazy

from ontology.models import OntologyTerm, OntologyService, OntologySnake, OntologyImportSource
from ontology.panel_app_ontology import update_gene_relations


class OntologyMatch:

    @dataclass
    class Score:
        name: str
        unit: float
        max: float
        note: str

        @property
        def score(self) -> float:
            return self.max * self.unit

        @property
        def percent(self) -> float:
            return self.unit * 100

        def as_json(self):
            return {
                "name": self.name,
                "max": self.max,
                "unit": self.unit,
                "score": self.score,
                "note": self.note
            }

    def __init__(self, term_id: str):
        self.term = OntologyTerm.get_or_stub(term_id)
        self.text_search: bool = False  # was this found from doing a text search on 3rd party search server or local database
        self.selected: bool = False  # has the user selected this term for whatever context this is
        self.direct_reference: bool = False  # was this term referenced by ID directly, e.g. text is "patient has MONDO:123456" and this is "MONDO:123456"
        self.gene_relationships: List[OntologySnake] = list()  # in what ways is this related to the gene in question (assuming there is a gene in question)

        self.scores: List[OntologyMatch.Score] = list()
        self.score = 0

    def __lt__(self, other):
        return self.score < other.score

    @staticmethod
    def _snake_as_json(snake: OntologySnake):
        via = None
        steps = snake.show_steps()
        if len(steps) > 1:
            via = steps[0].dest_term.id
        last_step = steps[0]
        return {
            "via": via,
            "extra": last_step.relation.extra,
            "relation": last_step.relation.relation,
            "source": last_step.relation.from_import.context
        }

    def as_json(self) -> Dict:
        data = {
            "id": self.term.id,
            "url": reverse("ontology_term", kwargs={"term": self.term.url_safe_id}),
            "title": self.term.name,
            "definition": self.term.definition,
            "scores": [score.as_json() for score in self.scores],
            "score": self.score,
            "gene_relationships": [OntologyMatch._snake_as_json(snake) for snake in self.gene_relationships]
        }
        if self.text_search:
            data["text_search"] = True
        if self.selected:
            data["selected"] = True
        if self.direct_reference:
            data["direct_reference"] = True

        return data


SKIP_TERMS = {"", "disease", "disorder", "the", "a", "an", "and", "or", "for", "the", "type"}
SUB_TYPE = re.compile("[0-9]+[a-z]?")


def tokenize_condition_text(text: str) -> Set[str]:
    if text is None:
        return set()
    text = text.replace(";", " ").replace("-", " ").lower()
    tokens = [token.strip() for token in text.split(" ")]
    tokens = [token for token in tokens if token not in SKIP_TERMS]
    return tokens


def pretty_set(s: Set[str]) -> str:
    tokens = list(s)
    tokens = sorted(tokens)
    tokens = [f'"{tok}"' for tok in tokens]
    text = ", ".join(tokens)
    return text


def normalize_condition_text(text: str):
    if text is None:
        return None
    text = text.lower()
    text = re.sub("[,;./]", " ", text)  # replace , ; . with spaces
    text = re.sub("[ ]{2,}", " ", text)  # replace multiple spaces with
    return text


class OntologyMatching:

    def __init__(self, search_term: Optional[str] = None, gene_symbol: Optional[str] = None):
        self.term_map: Dict[str, OntologyMatch] = dict()
        self.search_term = search_term
        self.search_terms = None
        if search_term:
            self.search_terms = set(tokenize_condition_text(normalize_condition_text(search_term))) - SKIP_TERMS

        self.gene_symbol = gene_symbol

    def find_or_create(self, term_id: str) -> OntologyMatch:
        mondo = self.term_map.get(term_id)
        if not mondo:
            mondo = OntologyMatch(term_id=term_id)
            self.term_map[term_id] = mondo

        return mondo

    def apply_scores(self):
        for match in self.term_map.values():
            self.apply_score(match)

    def apply_score(self, match: OntologyMatch):
        scores: List[OntologyMatch.Score] = list()
        regular_match = True
        if match.term.is_stub:
            # Should this be prioritised over direct reference?
            regular_match = False
            scores.append(OntologyMatch.Score(
                name="No copy of this term in our database",
                unit=1, max=-100
        ))
        if match.direct_reference:
            regular_match = False
            scores.append(OntologyMatch.Score(
                name="Directly referenced",
                unit=1, max=1000
        ))
        if regular_match:
            if search_terms := self.search_terms:
                match_terms = set(tokenize_condition_text(normalize_condition_text(match.term.name))) - SKIP_TERMS
                # TODO maybe check description for match terms?
                superfluous_words = match_terms.difference(search_terms)
                missing_words = search_terms.difference(match_terms)

                missing_word_ratio = 0 if not match_terms else (float(len(missing_words)) / float(len(search_terms)))
                superfluous_word_ratio = 0 if not match_terms else min(1, float(len(superfluous_words)) / float(len(match_terms)))

                # DIRECT_REFERENCE = OntologyScoreWeighting("Term is directly referenced", 1000)
                # GENE = OntologyScoreWeighting("Gene matching", 20)
                # WORD_MATCHING = OntologyScoreWeighting("Word matching", 40)
                # WORD_SUPERFLUOUS = OntologyScoreWeighting("Limited extra words bonus", 40)

                scores.append(OntologyMatch.Score(
                    name="Word Matching", max=40, unit=1 - missing_word_ratio,
                    note=f"Missing words {pretty_set(missing_words)}" if missing_words else "No missing words"
                ))
                scores.append(OntologyMatch.Score(
                    name="Limited extra words bonus", max=40, unit=1 - superfluous_word_ratio,
                    note=f"Superfluous words {pretty_set(superfluous_words)}" if superfluous_words else "No superfluous words"
                ))

            if gene_symbol := self.gene_symbol:
                # TODO work out if we want to apply different scores to different matches
                sources = set()
                for snake in match.gene_relationships:
                    relationship_source = snake.show_steps()[0].relation.from_import.import_source  # remember these steps are from gene to term, so the first step describes the relationship to the gene
                    sources.add(relationship_source)

                sources = list()
                if OntologyImportSource.PANEL_APP_AU in sources:
                    sources.append("PanelAPP AU")
                if OntologyImportSource.MONDO in sources:
                    sources.append("MONDO")
                if OntologyImportSource.HPO in sources:
                    sources.append("Entrez")

                scores.append(OntologyMatch.Score(
                    name="Gene relationship", max=20, unit=1 if sources else 0,
                    note=f"Term relates to gene according to {sources}" if sources else "No relationship between this term and gene in our database"
                ))
            else:
                scores.append(OntologyMatch.Score(
                    name="Gene relationship", max=20, unit=1,
                    note=f"Matching on gene not required at this level"
                ))
        match.scores = scores
        total = 0
        for score in scores:
            total += score.score
        match.score = total

    def populate_relationships(self):
        if gene_symbol := self.gene_symbol:
            update_gene_relations(gene_symbol)  # make sure panel app AU is up to date
            snakes = OntologySnake.terms_for_gene_symbol(gene_symbol=gene_symbol, desired_ontology=OntologyService.MONDO)  # always convert to MONDO for now
            seen_gene_relationships: Set[OntologyService] = set()  # if MONDO:123456 relates to OMIM:64454 via two different ways, we don't care
            # only care about the relationship to the gene
            for snake in snakes:
                gene_relation = snake.show_steps()[0].relation
                if gene_relation not in seen_gene_relationships:
                    seen_gene_relationships.add(gene_relation)
                    mondo_term = snake.leaf_term
                    mondo_meta = self.find_or_create(mondo_term.id)
                    mondo_meta.gene_relationships.append(snake)

    def select_term(self, term: str):
        self.find_or_create(term).selected = True

    def reference_term(self, term: str):
        self.find_or_create(term).direct_reference = True

    def searched_term(self, term: str):
        self.find_or_create(term).text_search = True

    def __iter__(self):
        values = sorted(list(self.term_map.values()), reverse=True)
        return iter(values)

    def as_json(self) -> Any:
        return [value.as_json() for value in self]

    @staticmethod
    def from_search(search_text: str, gene_symbol: Optional[str], selected: Optional[List[str]] = None, server_search: bool = True) -> 'OntologyMatching':
        ontology_matches = OntologyMatching(search_term=search_text, gene_symbol=gene_symbol)
        ontology_matches.populate_relationships()  # find all terms linked to the gene_symbol (if there is one)

        if selected:
            for select in selected:
                ontology_matches.select_term(select)

        ONTOLOGY_PATTERN = re.compile("(MONDO|OMIM|HP):[0-9]+")
        if ONTOLOGY_PATTERN.match(search_text):
            ontology_matches.find_or_create(search_text).direct_reference = True

        else:
            # Client search
            # this currently requires all words to be present
            search_terms = set(tokenize_condition_text(search_text))
            qs = OntologyTerm.objects.filter(ontology_service=OntologyService.MONDO)
            qs = qs.filter(reduce(operator.and_, [Q(name__icontains=term) for term in search_terms]))
            qs = qs.order_by(Length('name')).values_list("id", flat=True)
            result: OntologyTerm
            for result in qs[0:20]:
                ontology_matches.searched_term(result)

            # ALSO, do client side search
            if server_search:
                server_search_text = search_text
                if gene_symbol:
                    server_search_text = server_search_text + " " + gene_symbol
                try:
                    results = requests.get(f'https://api.monarchinitiative.org/api/search/entity/autocomplete/{server_search_text}', {
                        "prefix": "MONDO",
                        "rows": 6,
                        "minimal_tokenizer": "false",
                        "category": "disease"
                    }).json().get("docs")

                    for result in results:
                        o_id = result.get('id')
                        # result.get('label') gives the labe as it's known by the search server
                        ontology_matches.searched_term(o_id)
                except:
                    pass
                    # TODO communicate to the user couldn't search mondo text search

        ontology_matches.apply_scores()
        return ontology_matches
