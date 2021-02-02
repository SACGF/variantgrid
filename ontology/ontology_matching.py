import operator
import re
from dataclasses import dataclass
from functools import reduce
from typing import Dict, Optional, List, Any, Set, Iterable

import requests
from django.db.models import Q
from django.db.models.functions import Length
from django.urls import reverse
from lazy import lazy

from classification.regexes import db_ref_regexes, DbRegexes
from library.utils import empty_to_none
from ontology.models import OntologyTerm, OntologyService, OntologySnake, OntologyImportSource, OntologyTermRelation
from ontology.panel_app_ontology import update_gene_relations


class OntologyMatch:

    @dataclass
    class Score:
        name: str
        unit: float
        max: float
        note: Optional[str] = ""

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

    @lazy
    def is_leaf(self) -> Optional[bool]:
        if not self.term.is_stub and self.term.ontology_service == OntologyService.MONDO:
            return not OntologyTermRelation.children_of(self.term).exists()  # no children exist
        return None

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


SUFFIX_SKIP_TERMS = {"", "the", "an", "and", "or", "for", "the", "type", "group"}
PREFIX_SKIP_TERMS = SUFFIX_SKIP_TERMS.union({"a",})  # only exclude "A" from prefix, in case it says "type" A

SUB_TYPE = re.compile("^(.*?)(?: )((?:group|type)?(?: )?(?:[A-Z]|[0-9]+|[0-9]+[A-Z]|i|ii|iii|iv|v|vi|vii|viii|ix))$", re.IGNORECASE)
ROMAN = {
    "i": '1',
    "ii": '2',
    "iii": '3',
    "iv": '4',
    "v": '5',
    "vi": '6',
    "vii": '7',
    "viii": '8',
    "ix": '9'
}

class SearchText:

    @staticmethod
    def roman_to_arabic(numeral: str):
        return ROMAN.get(numeral, numeral)

    @staticmethod
    def tokenize_condition_text(text: str) -> Set[str]:
        if text is None:
            return set()
        text = text.replace(";", " ").replace("-", " ").lower()
        tokens = [token.strip() for token in text.split(" ")]
        return tokens

    def __init__(self, text: str):
        self.raw = text
        self.prefix = None
        self.prefix_terms: Set[str] = set()
        self.suffix = None
        self.suffix_terms: Set[str] = set()

        normal_text = normalize_condition_text(text)
        if sub_type_match := SUB_TYPE.match(normal_text):
            self.prefix = sub_type_match.group(1).strip()
            self.suffix = sub_type_match.group(2).strip()

            self.prefix_terms = set(SearchText.tokenize_condition_text(self.prefix)) - PREFIX_SKIP_TERMS
            self.suffix_terms = set(SearchText.roman_to_arabic(term) for term in SearchText.tokenize_condition_text(self.suffix)) - SUFFIX_SKIP_TERMS
        else:
            self.prefix = normal_text
            self.prefix_terms = set(SearchText.tokenize_condition_text(self.prefix)) - PREFIX_SKIP_TERMS

    @property
    def prefix_terms_display(self) -> str:
        return ", ".join(sorted(list(self.prefix_terms)))

    @property
    def suffix_terms_display(self):
        # TODO make these django filters
        return ", ".join(sorted(list(self.suffix_terms)))

    @property
    def all_terms(self) -> Set[str]:
        return self.prefix_terms.union(self.suffix_terms)


def pretty_set(s: Iterable[str]) -> str:
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
    text = text.strip()
    return text


class OntologyMatching:

    def __init__(self, search_term: Optional[str] = None, gene_symbol: Optional[str] = None):
        self.term_map: Dict[str, OntologyMatch] = dict()
        self.search_text: Optional[SearchText] = None
        self.sub_type = None
        if search_term:
            self.search_text = SearchText(search_term)
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
        if match.direct_reference:
            regular_match = False
            scores.append(OntologyMatch.Score(
                name="Directly referenced",
                unit=1, max=1000, note=""
            ))
        if match.term.is_stub:
            # Should this be prioritised over direct reference?
            regular_match = False
            scores.append(OntologyMatch.Score(
                name="No copy of this term in our database",
                unit=1, max=-100, note=""
            ))
        if match.term.is_obsolete:
            regular_match = False
            scores.append(OntologyMatch.Score(
                name="This term is marked as obsolete",
                unit=1, max=-100, note=""
            ))

        if regular_match:
            match_text = SearchText(match.term.name)
            if search_text := self.search_text and self.search_text.raw:

                superfluous_words = set()
                missing_words = set()
                missing_word_ratio = 0
                superfluous_word_ratio = 0

                search_text_terms: Set[str]
                match_text_terms: Set[str]
                if self.search_text.suffix_terms or not self.gene_symbol:
                    search_text_terms = self.search_text.all_terms
                    match_text_terms = match_text.all_terms
                else:
                    # match has a suffix and search term doesn't (but does have a gene symbol)
                    # so don't penalise for the suffix
                    search_text_terms = self.search_text.prefix_terms
                    match_text_terms = match_text.prefix_terms

                missing_words = search_text_terms.difference(match_text_terms)
                missing_word_ratio = float(len(missing_words)) / float(len(search_text_terms))

                superfluous_words = match_text_terms.difference(search_text_terms)
                superfluous_word_ratio = float(len(superfluous_words)) / float(len(match_text_terms))

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
                source_codes = set()
                for snake in match.gene_relationships:
                    relationship_source = snake.show_steps()[0].relation.from_import.import_source  # remember these steps are from gene to term, so the first step describes the relationship to the gene
                    source_codes.add(relationship_source)

                sources = list()
                # these are the only sources we care about?
                if OntologyImportSource.PANEL_APP_AU in source_codes:
                    sources.append("PanelAPP AU")
                if OntologyImportSource.MONDO in source_codes:
                    sources.append("MONDO")
                if OntologyImportSource.HPO in source_codes:
                    sources.append("NCBI")

                if not sources:
                    scores.append(OntologyMatch.Score(
                        name="Gene relationship", max=20, unit=0,
                        note="No relationship between this term and gene in our database"
                    ))
                else:
                    scores.append(OntologyMatch.Score(
                        name="Gene relationship", max=20, unit=1,
                        note=f"Has relationships via {pretty_set(sources)}"
                    ))
                    if not match.is_leaf:
                        scores.append(OntologyMatch.Score(
                            name="Gene relationship - not leaf leaf term", max=-2, unit=1,
                            note="Term has children, could be more specific"
                        ))
            else:
                scores.append(OntologyMatch.Score(
                    name="Gene relationship", max=20, unit=1,
                    note=f"Matching on gene not required at this level"
                ))

            injected_suffix = match_text.suffix_terms and not self.search_text.suffix_terms
            if injected_suffix:
                scores.append(OntologyMatch.Score(
                    name="Gene relationship - Extra suffix", max=-1, unit=1,
                    note=f"Has extra suffix of '{match_text.suffix}'"
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
            for snake in snakes:
                gene_relation = snake.show_steps()[0].relation
                mondo_term = snake.leaf_term
                mondo_meta = self.find_or_create(mondo_term.id)
                mondo_meta.gene_relationships.append(snake)

    def select_term(self, term: str):
        self.find_or_create(term).selected = True

    def reference_term(self, term: str):
        self.find_or_create(term).direct_reference = True

    def searched_term(self, term: str):
        self.find_or_create(term).text_search = True

    def top_terms(self) -> List[OntologyMatch]:
        top_values: List[OntologyMatch] = list()
        for match in self:
            if len(top_values) == 0:
                top_values.append(match)
            elif top_values[0].score == match.score:
                top_values.append(match)
            else:
                return top_values
        return top_values

    def __iter__(self):
        values = sorted(list(self.term_map.values()), reverse=True)
        return iter(values)

    def as_json(self) -> Any:
        return [value.as_json() for value in self]

    @staticmethod
    def from_search(search_text: str, gene_symbol: Optional[str], selected: Optional[List[str]] = None, server_search: bool = True) -> 'OntologyMatching':
        search_text = empty_to_none(search_text)
        ontology_matches = OntologyMatching(search_term=search_text, gene_symbol=gene_symbol)
        ontology_matches.populate_relationships()  # find all terms linked to the gene_symbol (if there is one)

        if selected:
            for select in selected:
                ontology_matches.select_term(select)

        if search_text:
            matches = db_ref_regexes.search(search_text, default_regex=DbRegexes.OMIM)
            has_match = not not matches
            for match in matches:
                ontology_matches.find_or_create(match.id_fixed).direct_reference = True

            if not has_match:
                # only search if we didn't find an OMIM/MONDO term linked directly

                # Client search this currently requires all words to be present
                search_terms = set(SearchText.tokenize_condition_text(search_text))
                qs = OntologyTerm.objects.filter(ontology_service=OntologyService.MONDO)
                qs = qs.filter(reduce(operator.and_, [Q(name__icontains=term) for term in search_terms]))
                qs = qs.order_by(Length('name')).values_list("id", flat=True)
                result: OntologyTerm
                for result in qs[0:20]:
                    ontology_matches.searched_term(result)

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
