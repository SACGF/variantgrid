import operator
import re
from dataclasses import dataclass
from functools import reduce
from typing import Dict, Optional, List, Any, Set, Iterable

import requests
from django.db.models import Q
from django.urls import reverse
from lazy import lazy

from classification.regexes import db_ref_regexes
from library.log_utils import report_message
from library.utils import empty_to_none
from ontology.models import OntologyTerm, OntologyService, OntologySnake
from ontology.panel_app_ontology import update_gene_relations


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
SUFFIX_SKIP_TERMS = {"", "the", "an", "and", "or", "for", "the", "type", "group"}
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
    "ix": '9'
}

@dataclass
class MatchInfo:
    alias_index: Optional[int] = None


class SearchText:

    @staticmethod
    def roman_to_arabic(numeral: str):
        return ROMAN.get(numeral, numeral)

    @staticmethod
    def tokenize_condition_text(text: str, deplural=False) -> Set[str]:
        if text is None:
            return set()
        text = text.replace(";", " ").replace("-", " ").lower()
        tokens = [token.strip() for token in text.split(" ")]
        if deplural:
            new_tokens = list()
            for token in tokens:
                if len(token) >= 5 and token.endswith('s'):  # make sure 5 or more characters long so not acronymn
                    token = token[0:-1]
                new_tokens.append(token)
            tokens = new_tokens
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

            self.prefix_terms = set(SearchText.tokenize_condition_text(self.prefix, deplural=True)) - PREFIX_SKIP_TERMS
            self.suffix_terms = set(SearchText.roman_to_arabic(term) for term in SearchText.tokenize_condition_text(self.suffix)) - SUFFIX_SKIP_TERMS
        else:
            self.prefix = normal_text
            self.prefix_terms = set(SearchText.tokenize_condition_text(self.prefix, deplural=True)) - PREFIX_SKIP_TERMS

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
        return self.prefix_terms == other.prefix_terms and self.suffix_terms == other.suffix_terms


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


@dataclass
class TopLevelSuggestionDetails:
    terms: Optional[List[OntologyTerm]]
    ids_in_text: bool = False
    ids_in_text_warning: str = None
    checked_local_terms: int = 0
    checked_server_terms: int = 0
    alias_index: Optional[int] = None

    def is_auto_assignable(self):
        if terms := self.terms:
            if len(terms) != 1:
                return False
            if self.ids_in_text_warning:
                return False
            if self.alias_index is not None and self.alias_index > 1:
                return False
            if not self.is_leaf:
                return False
            if self.is_obsolete:
                return False
            return True
        return False

    def __bool__(self):
        return bool(self.terms)

    @property
    def is_obsolete(self):
        if terms := self.terms:
            for term in terms:
                if not term.is_obsolete:
                    return True
        return False

    def is_leaf(self):
        if terms := self.terms:
            for term in terms:
                if not term.is_leaf:
                    return False
            return True
        return None


class OntologyMatching:

    @staticmethod
    def top_level_suggestion(text: str, fallback_to_online: bool = True) -> TopLevelSuggestionDetails:
        if suggestion := OntologyMatching.embedded_ids_check(text):
            return suggestion
        return OntologyMatching.search_suggestion(text, fallback_to_online=fallback_to_online)

    @staticmethod
    def embedded_ids_check(text: str) -> Optional[TopLevelSuggestionDetails]:

        db_matches = db_ref_regexes.search(text)
        detected_any_ids = bool(db_matches)  # see if we found any prefix suffix, if we do,
        db_matches = [match for match in db_matches if match.db in ["OMIM", "HP", "MONDO"]]

        found_terms = list()
        for match in db_matches:
            found_terms.append(OntologyTerm.get_or_stub(match.id_fixed))

        found_stray_omim = False

        # fall back to looking for stray OMIM terms if we haven't found any ids e.g. PMID:123456 should stop this code
        if not detected_any_ids:
            stray_omim_matches = OPRPHAN_OMIM_TERMS.findall(text)
            stray_omim_matches = [term for term in stray_omim_matches if len(term) == 6]
            if stray_omim_matches:
                for omim_index in stray_omim_matches:
                    omim = OntologyTerm.get_or_stub(f"OMIM:{omim_index}")
                    found_terms.append(omim)
                    found_stray_omim = True
        if found_terms:
            warning = None
            if found_stray_omim:
                warning = "Detected OMIM terms from raw numbers without prefix"

            return TopLevelSuggestionDetails(terms=found_terms, ids_in_text=True, ids_in_text_warning=warning)

        return None

    @staticmethod
    def search_suggestion(text: str, fallback_to_online: bool = True) -> TopLevelSuggestionDetails:
        match_text = SearchText(text)
        q = list()
        # TODO, can we leverage phenotype matching?
        for term in match_text.prefix_terms:
            if len(term) > 2:
                # TODO evaluate if it was worth it comparing aliases
                q.append(Q(name__icontains=term) | Q(aliases__icontains=term))
        # don't bother with searching for suffix, just find them all and see how we go with the matching
        local_term_count = 0
        if q:

            for term in OntologyTerm.objects.filter(ontology_service__in={OntologyService.MONDO, OntologyService.OMIM}).filter(reduce(operator.and_, q)).order_by('ontology_service')[0:100]:
                local_term_count += 1
                if match_info := match_text.matches(term):
                    return TopLevelSuggestionDetails(
                        terms=[term],
                        checked_local_terms=local_term_count,
                        alias_index=match_info.alias_index
                    )

        server_term_count = 0
        if fallback_to_online:
            try:
                # TODO ensure "text" is safe, it should already be normalised
                results = requests.get(
                    f'https://api.monarchinitiative.org/api/search/entity/autocomplete/{text}', {
                        "prefix": "MONDO",
                        "rows": 10,
                        "minimal_tokenizer": "false",
                        "category": "disease"
                    }).json().get("docs")

                for result in results:
                    server_term_count += 1
                    o_id = result.get('id')
                    # result.get('label') gives the label as it's known by the search server
                    term = OntologyTerm.get_or_stub(o_id)
                    if match_info := match_text.matches(term):
                        return TopLevelSuggestionDetails(
                            terms=[term],
                            checked_local_terms=local_term_count,
                            checked_server_terms=server_term_count,
                            alias_index=match_info.alias_index
                        )
            except:
                print("Error searching server")

        return TopLevelSuggestionDetails(
            terms=None,
            checked_local_terms=local_term_count,
            checked_server_terms=server_term_count,
        )

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

    """
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
                injected_suffix = match_text.suffix_terms and not self.search_text.suffix_terms
                if injected_suffix:
                    scores.append(OntologyMatch.Score(
                        name="Extra suffix", max=-1, unit=1,
                        note=f"Has extra suffix of '{match_text.suffix}'"
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

        match.scores = scores
        total = 0
        for score in scores:
            total += score.score
        match.score = total
        
    def only_top_term(self) -> Optional[OntologyMatch]:
        tops = self.top_terms()
        if len(tops) == 1:
            return tops[0]
        else:
            return None

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
    """

    def populate_relationships(self, server_search=True):
        if gene_symbol := self.gene_symbol:
            try:
                OntologyTerm.get_gene_symbol(gene_symbol)
            except ValueError:
                report_message(message=f"Could not resolve {gene_symbol} to HGNC OntologyTerm - can't do gene specific resolutions", level='warning')
                return

            if server_search:
                update_gene_relations(gene_symbol)  # make sure panel app AU is up to date
            snakes = OntologySnake.terms_for_gene_symbol(gene_symbol=gene_symbol, desired_ontology=OntologyService.MONDO)  # always convert to MONDO for now
            for snake in snakes:
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
        return [value.as_json() for value in self]

    @staticmethod
    def from_search(search_text: str, gene_symbol: Optional[str], selected: Optional[List[str]] = None) -> 'OntologyMatching':
        search_text = empty_to_none(search_text)
        ontology_matches = OntologyMatching(search_term=search_text, gene_symbol=gene_symbol)
        ontology_matches.populate_relationships(server_search=True)  # find all terms linked to the gene_symbol (if there is one)

        if selected:
            for select in selected:
                ontology_matches.select_term(select)

        if search_text:
            matches = db_ref_regexes.search(search_text)
            detected_any_ids = not not matches
            detected_ontology_id = False
            matches = [match for match in matches if match.db in ["OMIM", "HP", "MONDO"]]

            for match in matches:
                detected_ontology_id = True
                ontology_matches.find_or_create(match.id_fixed).direct_reference = True

            # fall back to looking for stray OMIM terms if we haven't found any ids e.g. PMID:123456 should stop this code
            if not detected_any_ids:
                stray_omim_matches = OPRPHAN_OMIM_TERMS.findall(search_text)
                stray_omim_matches = [term for term in stray_omim_matches if len(term) == 6]
                if stray_omim_matches:
                    detected_ontology_id = True
                    for omim_index in stray_omim_matches:
                        ontology_matches.find_or_create(f"OMIM:{omim_index}").direct_reference = True

            if not detected_ontology_id:
                """
                # technically this bit isn't a server search, but need to make the boolean more flexible
                search_terms = set(SearchText.tokenize_condition_text(search_text))
                qs = OntologyTerm.objects.filter(ontology_service=OntologyService.MONDO)
                qs = qs.filter(reduce(operator.and_, [Q(name__icontains=term) for term in search_terms]))
                qs = qs.order_by(Length('name')).values_list("id", flat=True)
                result: OntologyTerm
                for result in qs[0:20]:
                    ontology_matches.searched_term(result)
                """
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
