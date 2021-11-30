import re
from dataclasses import dataclass
from typing import Dict, Optional, List, Any, Set, Iterable, TypedDict

import requests
from django.urls import reverse
from lazy import lazy

from annotation.regexes import db_ref_regexes
from library.log_utils import report_message
from library.utils import empty_to_none
from ontology.models import OntologyTerm, OntologyService, OntologySnake, OntologyImportSource


class OntologySnakeJson(TypedDict):
    via: str
    extra: Optional[Dict[Any, Any]]
    relation: str
    source: str


class OntologyMatch:

    def __init__(self, term_id: str):
        self.term = OntologyTerm.get_or_stub(term_id)
        self.selected: bool = False  # has the user selected this term for whatever context this is
        self.direct_reference: bool = False  # was this term referenced by ID directly, e.g. text is "patient has MONDO:123456" and this is "MONDO:123456"
        self.gene_relationships: List[OntologySnake] = list()  # in what ways is this related to the gene in question (assuming there is a gene in question)
        self.search_engine_score: Optional[int] = None  # was this found from doing a text search on 3rd party search, and if so, what's its ranking

    @property
    def text_search(self):
        return self.search_engine_score is not None

    @lazy
    def _score(self):
        score = 0
        if self.selected:
            score += 10000
        if self.direct_reference:
            score += 1000
        if bool(self.gene_relationships):
            score += 100
        if self.search_engine_score:
            score += self.search_engine_score
        return score

    @property
    def is_leaf(self) -> Optional[bool]:
        return self.term.is_leaf

    def __lt__(self, other: 'OntologyMatch'):
        return self._score < other._score

    @staticmethod
    def _snake_as_json(snake: OntologySnake) -> OntologySnakeJson:
        via = None
        steps = snake.show_steps()
        if len(steps) > 1:
            via = steps[0].dest_term.id
        last_step = steps[0]
        return OntologySnakeJson(
            via=via,
            extra=last_step.relation.extra,
            relation=last_step.relation.relation,
            source=last_step.relation.from_import.context
        )

    def as_json(self) -> Dict:

        data = {
            "id": self.term.id,
            "url": reverse("ontology_term", kwargs={"term": self.term.url_safe_id}),
            "title": self.term.name,
            "definition": self.term.definition,
            "gene_relationships": [OntologyMatch._snake_as_json(snake) for snake in self.gene_relationships],
            "obsolete": self.term.is_obsolete
        }
        if self.text_search:
            data["text_search"] = True
        if self.selected:
            data["selected"] = True
        if self.direct_reference:
            data["direct_reference"] = True

        return data


OPRPHAN_OMIM_TERMS = re.compile("[0-9]{6,}")
SUFFIX_SKIP_TERMS = {"", "the", "an", "and", "or", "for", "the", "type", "group", "with"}
PREFIX_SKIP_TERMS = SUFFIX_SKIP_TERMS.union({"a", })  # only exclude "A" from prefix, in case it says "type" A

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
    "ix": '9',
    # 'x': '10', matching on x can be dangerous, has a lot of other meanings
    'xi': '11',
    'xii': '12',
    'xiii': '13',
    'xiv': '14',
    'xv': '15',
    'xvi': '16',
    'xvii': '17',
    'xviii': '18',
    'ixx': '19'
}


@dataclass
class MatchInfo:
    alias_index: Optional[int] = None
    """ What alias (if any) did we use to match, None if matched on name or embedded ID """


class SearchText:  # TODO shold be renamed ConditionSearchText
    """
    Text broken up into tokens, optionally de-pluralised, de-roman numeraled
    split into prefix and suffix (if we detect a splitter word like 'type')
    see the SUB_TYPE regex for details
    """

    @staticmethod
    def roman_to_arabic(numeral: str):
        return ROMAN.get(numeral, numeral)

    @staticmethod
    def tokenize_condition_text(text: str, deplural=False, deroman=False) -> Set[str]:
        if text is None:
            return set()
        text = text.replace(";", " ").replace("-", " ").lower()
        tokens = [token.strip() for token in text.split(" ")]
        if deplural:
            new_tokens = list()
            for token in tokens:
                if len(token) >= 5 and token.endswith('s'):  # make sure 5 or more characters long so not acronym
                    token = token[0:-1]
                new_tokens.append(token)
            tokens = new_tokens
        if deroman:
            tokens = [SearchText.roman_to_arabic(token) for token in tokens]
        return set(tokens)

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

            self.prefix_terms = set(SearchText.tokenize_condition_text(self.prefix, deroman=True, deplural=True)) - PREFIX_SKIP_TERMS
            self.suffix_terms = set(SearchText.tokenize_condition_text(self.suffix, deroman=True)) - SUFFIX_SKIP_TERMS
        else:
            self.prefix = normal_text
            self.prefix_terms = set(SearchText.tokenize_condition_text(self.prefix, deroman=True, deplural=True)) - PREFIX_SKIP_TERMS

    @property
    def prefix_terms_display(self) -> str:
        return ", ".join(sorted(list(self.prefix_terms)))

    @property
    def suffix_terms_display(self):
        return ", ".join(sorted(list(self.suffix_terms)))

    @property
    def all_terms(self) -> Set[str]:
        return self.prefix_terms.union(self.suffix_terms)

    def matches(self, term: OntologyTerm) -> Optional[MatchInfo]:
        if name := term.name:
            if SearchText(name).effective_equals(self):
                return MatchInfo(None)
        if aliases := term.aliases:
            for index, alias in enumerate(aliases):
                if SearchText(alias).effective_equals(self):
                    return MatchInfo(index)
        return None

    def effective_equals(self, other: 'SearchText') -> bool:
        # TODO handle a little bit off by 1 letter matching
        return self.prefix_terms == other.prefix_terms and self.suffix_terms == other.suffix_terms or self.all_terms == other.all_terms


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
    text = re.sub("[,;./?()]", " ", text)  # replace , ; . with spaces
    text = re.sub("[ ]{2,}", " ", text)  # replace multiple spaces with
    text = text.strip()
    return text


class OntologyMatching:
    """
    Used to build up a big list of suggestions, typically on the MondoPicker
    """

    def __init__(self, search_term: Optional[str] = None, gene_symbol: Optional[str] = None):
        self.term_map: Dict[str, OntologyMatch] = dict()
        self.errors: List[str] = list()
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

    def add_error(self, error: str):
        self.errors.append(error)

    def populate_relationships(self):
        """
        Given a gene symbol, provide all terms that have a relationship to that gene symbol
        """
        if gene_symbol := self.gene_symbol:
            try:
                OntologyTerm.get_gene_symbol(gene_symbol)
            except ValueError:
                report_message(message=f"Could not resolve {gene_symbol} to HGNC OntologyTerm - can't do gene specific resolutions", level='warning')
                return

            snakes = OntologySnake.terms_for_gene_symbol(gene_symbol=gene_symbol, desired_ontology=OntologyService.MONDO)  # always convert to MONDO for now
            had_panel_app = False
            for snake in snakes:
                if snake.show_steps()[0].relation.from_import.import_source == OntologyImportSource.PANEL_APP_AU:
                    if had_panel_app:
                        continue
                    had_panel_app = True

                mondo_term = snake.leaf_term
                mondo_meta = self.find_or_create(mondo_term.id)
                mondo_meta.gene_relationships.append(snake)  # assign the snake to the term

    def select_term(self, term: str):
        self.find_or_create(term).selected = True

    def reference_term(self, term: str):
        self.find_or_create(term).direct_reference = True

    def searched_term(self, term: str, score: int):
        self.find_or_create(term).search_engine_score = score

    def __iter__(self):
        values = sorted(list(self.term_map.values()), reverse=True)
        return iter(values)

    def as_json(self) -> Any:
        return {"errors": self.errors, "terms": [value.as_json() for value in self]}

    @staticmethod
    def from_search(search_text: str, gene_symbol: Optional[str], selected: Optional[List[str]] = None) -> 'OntologyMatching':
        search_text = empty_to_none(search_text)
        ontology_matches = OntologyMatching(search_term=search_text, gene_symbol=gene_symbol)
        ontology_matches.populate_relationships()  # find all terms linked to the gene_symbol (if there is one)

        if selected:
            for select in selected:
                try:
                    ontology_matches.select_term(select)
                except ValueError:
                    ontology_matches.add_error(f"\"{select}\" is not a valid ontology term.")

        if search_text:
            # WARNING, this code is duplicated in condition_text_matching.embedded_ids_check
            matches = db_ref_regexes.search(search_text)
            detected_any_ids = not not matches
            detected_ontology_id = False
            matches = [match for match in matches if match.db in OntologyService.CONDITION_ONTOLOGIES]

            for match in matches:
                detected_ontology_id = True
                ontology_matches.find_or_create(match.id_fixed).direct_reference = True

            # fall back to looking for stray OMIM terms if we haven't found any ids e.g. PMID:123456 should stop this code
            if not detected_any_ids and ':' not in search_text:
                stray_omim_matches = OPRPHAN_OMIM_TERMS.findall(search_text)
                stray_omim_matches = [term for term in stray_omim_matches if len(term) == 6]
                if stray_omim_matches:
                    detected_ontology_id = True
                    for omim_index in stray_omim_matches:
                        ontology_matches.find_or_create(f"OMIM:{omim_index}").direct_reference = True

            if not detected_ontology_id:
                # the actual server search
                server_search_text = search_text
                if gene_symbol:
                    server_search_text = server_search_text + " " + gene_symbol
                try:
                    row_count = 6
                    results = requests.get(f'https://api.monarchinitiative.org/api/search/entity/autocomplete/{server_search_text}', {
                        "prefix": "MONDO",
                        "rows": row_count,
                        "minimal_tokenizer": "false",
                        "category": "disease"
                    }).json().get("docs")

                    for index, result in enumerate(results):
                        o_id = result.get('id')
                        # result.get('label') gives the label as it's known by the search server
                        ontology_matches.searched_term(o_id, row_count - index)
                except:
                    pass
                    # TODO communicate to the user couldn't search mondo text search

        # ontology_matches.apply_scores()
        return ontology_matches
