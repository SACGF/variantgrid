import operator
import re
from dataclasses import dataclass
from functools import reduce
from typing import Dict, Optional, List, Any, Set

import requests
from django.db.models import Q
from django.db.models.functions import Length
from django.urls import reverse

from ontology.models import OntologyTerm, OntologyTermRelation, OntologyService, OntologySnake, OntologyRelation
from ontology.panel_app_ontology import update_gene_relations


class OntologyContext:

    @property
    def name(self):
        raise NotImplementedError("Ontology context must implement name")

    def merge(self, merge_me):
        pass

    def as_json(self) -> Dict:
        raise NotImplementedError("Ontology context must implement as_json")


class MetaKeys:
    SEARCHED = "searched"
    REFERENCED = "referenced"
    SELECTED = "selected"

    MONDO_GENE_LINK = "Monarch"
    PANELAPP_GENE_LINK = "PanelApp AU"
    CHILD_GENE_LINK = "Child"


class OntologyContextSearched(OntologyContext):
    """
    Indicates that this was found from doing a text search
    e.g. search for "rett syndrome" and found "MONDO:0010726" based on text match.
    """

    @property
    def name(self):
        return MetaKeys.SEARCHED

    # TODO might include search score in future
    def as_json(self) -> Any:
        return True


class OntologyContextDirectReference(OntologyContext):
    """
    Indicates that the search text actually made explicit reference to this very term
    e.g. "Patient suffers from MONDO:00000343"
    """

    @property
    def name(self):
        return MetaKeys.REFERENCED

    def as_json(self) -> Any:
        return True


class OntologyContextSelected(OntologyContext):
    """
    Indicates that this value has already been selected for the text
    """

    @property
    def name(self):
        return MetaKeys.SELECTED

    def as_json(self) -> Any:
        return True


class OntologyContextMonarchLink(OntologyContext):
    """
    Indicates that the mondo file determined a relationship between the given gene symbol
    and this ontology term
    """

    def __init__(self, gene_symbol: str, relation: str):
        self.gene_symbol = gene_symbol
        self.relation = relation

    @property
    def name(self):
        return MetaKeys.MONDO_GENE_LINK

    def as_json(self) -> Dict:
        return {
            "gene_symbol": self.gene_symbol,
            "relation": self.relation
        }


class OntologyContextPanelApp(OntologyContext):
    """
    Indicates that PanelApp has a known relationship between the given gene symbol
    and this ontology term
    """

    def __init__(self, gene_symbol: str, omim_id: int, phenotype_row: str, evidence: List[str]):
        self.gene_symbol = gene_symbol
        self.omim_id = omim_id
        self.phenotype_row = phenotype_row
        self.evidence = evidence

    def merge(self, merge_me):
        self.evidence = self.evidence + merge_me.evidence

    @property
    def name(self):
        return MetaKeys.PANELAPP_GENE_LINK

    def as_json(self):
        sorted_evidence = list(set(self.evidence))
        sorted_evidence.sort()
        return {
            "gene_symbol": self.gene_symbol,
            "omim_id": self.omim_id,
            "phenotype_row": self.phenotype_row,
            "evidence": sorted_evidence,
        }


class OntologContextChildGeneRelationship(OntologyContext):
    """
    Indicates that a child of this term has relationship to the gene symbol
    """

    def __init__(self, child_term_id: str):
        self.children = [child_term_id]

    def merge(self, merge_me):
        self.children += merge_me.children

    @property
    def name(self):
        return MetaKeys.CHILD_GENE_LINK

    def as_json(self):
        return {
            "matching_children": self.children,
        }

@dataclass
class OntologyScoreWeighting:
    name: str
    weighting: float


class OntologyMatchScoreWeightings:
    DIRECT_REFERENCE = OntologyScoreWeighting("Term is directly referenced", 1000)
    GENE = OntologyScoreWeighting("Gene matching", 20)
    WORD_MATCHING = OntologyScoreWeighting("Word matching", 40)
    WORD_SUPERFLUOUS = OntologyScoreWeighting("Limited extra words bonus", 40)


@dataclass
class OntologyMatchScorePart:
    weighting: OntologyScoreWeighting
    unit: float
    note: str

    @property
    def score(self) -> float:
        return self.weighting.weighting * self.unit

    @property
    def name(self) -> str:
        return self.weighting.name

    @property
    def max(self) -> float:
        return self.weighting.weighting

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


class OntologyMeta:

    def __init__(self, term_id: str):
        self.term = OntologyTerm.get_or_stub(term_id)
        self.scores: List[OntologyMatchScorePart] = list()
        self.contexts: Dict[str, OntologyContext] = dict()

    @property
    def score(self) -> float:
        score = 0
        for score_part in self.scores:
            score += score_part.score
        return score

    def __lt__(self, other):
        return self.score < other.score

    def add_score(self, score: OntologyMatchScorePart):
        self.scores.append(score)

    def add_context(self, context: OntologyContext):
        if existing := self.contexts.get(context.name):
            existing.merge(context)
        else:
            self.contexts[context.name] = context

    def as_json(self) -> Dict:
        sorted_names = list(self.contexts.keys())
        sorted_names.sort()
        contexts = {key: self.contexts[key].as_json() for key in sorted_names}
        scores = [score.as_json() for score in self.scores]

        data = {
            "id": self.term.id,
            "url": reverse("ontology_term", kwargs={"term": self.term.url_safe_id}),
            "title": self.term.name,
            "definition": self.term.definition,
            "contexts": contexts,
            "scores": scores,
            "score": self.score
        }

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
        self.term_map: Dict[str, OntologyMeta] = dict()
        self.search_term = search_term
        self.search_terms = None
        if search_term:
            self.search_terms = set(tokenize_condition_text(normalize_condition_text(search_term))) - SKIP_TERMS

        self.gene_symbol = gene_symbol

    def find_or_create(self, term_id: str) -> OntologyMeta:
        mondo = self.term_map.get(term_id)
        if not mondo:
            mondo = OntologyMeta(term_id=term_id)
            self.term_map[term_id] = mondo

            if self.gene_symbol:
                pass
            else:
                mondo.add_score(OntologyMatchScorePart(
                    OntologyMatchScoreWeightings.GENE,
                    unit=1,
                    note="Not matching at gene level, gene association not required/relevant"
                ))

        return mondo

    def apply_word_match_score(self):
        for onto in self.term_map.values():
            self.add_word_match_score(onto)

    def add_word_match_score(self, mondo: OntologyMeta):
        if search_terms := self.search_terms:
            if mondo.term.is_stub:
                mondo.add_score(OntologyMatchScorePart(
                    OntologyScoreWeighting("No copy of this term in our database", -100),
                    unit=-1
                ))
            else:
                match_terms = set(tokenize_condition_text(normalize_condition_text(mondo.term.name))) - SKIP_TERMS
                superfluous_words = match_terms.difference(search_terms)
                missing_words = search_terms.difference(match_terms)

                missing_word_ratio = 0 if not match_terms else (float(len(missing_words)) / float(len(search_terms)))
                superfluous_word_ratio = 0 if not match_terms else min(1, float(len(superfluous_words)) / float(len(match_terms)))

                mondo.add_score(OntologyMatchScorePart(
                    OntologyMatchScoreWeightings.WORD_MATCHING,
                    unit=1 - missing_word_ratio,
                    note=f"Missing words {pretty_set(missing_words)}" if missing_words else "No missing words"
                ))

                mondo.add_score(OntologyMatchScorePart(
                    OntologyMatchScoreWeightings.WORD_SUPERFLUOUS,
                    unit=1 - superfluous_word_ratio,
                    note=f"Superfluous words {pretty_set(superfluous_words)}" if superfluous_words else "No superfluous words",
                ))

                # if superfluous_words and not self.gene_symbol:
                #     for word in superfluous_words:
                #         if SUB_TYPE.match(word):
                #             mondo.add_score(OntologyMatchScorePart(
                #                 "Superfluous words includes subtype \"{word}\"",
                #                 20.0
                #             ))

    def populate_relationships(self):
        scores: Dict[str, OntologyMeta] = dict()
        if gene_symbol := self.gene_symbol:
            update_gene_relations(gene_symbol=gene_symbol)
            snakes = OntologySnake.terms_for_gene_symbol(gene_symbol=gene_symbol, desired_ontology=OntologyService.MONDO)
            for snake in snakes:
                # always convert to MONDO for now
                mondo_term = snake.leaf_term
                mondo_meta = self.find_or_create(mondo_term.id)
                scores[mondo_term.id] = mondo_meta
                gene_relation = snake.paths[0]
                gene_relationships_via = list()
                if gene_relation.relation == OntologyRelation.PANEL_APP_AU:
                    mondo_meta.add_context(OntologyContextPanelApp(
                        gene_symbol=gene_symbol,
                        omim_id=0, # TODO populate the OMIM ID
                        phenotype_row=gene_relation.extra.get("phenotype_row"),
                        evidence=gene_relation.extra.get("evidence")
                    ))
                else:
                    mondo_meta.add_context(OntologyContextMonarchLink(
                        gene_symbol=gene_symbol,
                        relation=gene_relation.relation
                    ))

                for parent_term in OntologyTermRelation.parents_of(mondo_term):
                    parent_meta = self.find_or_create(parent_term.id)
                    scores[parent_term.id] = parent_meta
                    if not MetaKeys.MONDO_GENE_LINK in parent_meta.contexts and \
                        not MetaKeys.PANELAPP_GENE_LINK in parent_meta.contexts:
                        parent_meta.add_context(OntologContextChildGeneRelationship(
                            child_term_id=mondo_term.id
                        ))

        for meta in scores.values():
            gene_relationships_via = set()
            if MetaKeys.MONDO_GENE_LINK in meta.contexts:
                gene_relationships_via.add("MONDO data")
            if MetaKeys.PANELAPP_GENE_LINK in meta.contexts:
                gene_relationships_via.add("PanelApp AU")
            if gene_relationships_via:
                gene_relationships_via_str = ", ".join(gene_relationships_via)
                meta.add_score(OntologyMatchScorePart(
                    OntologyMatchScoreWeightings.GENE,
                    unit=1,
                    note=f"Established gene relationship (through {gene_relationships_via_str})"
                ))
            elif MetaKeys.CHILD_GENE_LINK in meta.contexts:
                # TODO should this be removed now that we establish parent MONDO relationships
                meta.add_score(OntologyMatchScorePart(
                    OntologyMatchScoreWeightings.GENE,
                    unit=0.5,
                    note=f"Parent of term that has gene relationship"
                ))
            else:
                meta.add_score(OntologyMatchScorePart(
                    OntologyMatchScoreWeightings.GENE,
                    unit=0,
                    note=f"No relationship to gene symbol in database"
                ))

    def select_term(self, term: str):
        self.find_or_create(term).add_context(OntologyContextSelected())

    def reference_term(self, term: str):
        onto_meta = self.find_or_create(term)
        onto_meta.add_score(OntologyMatchScorePart(
            OntologyMatchScoreWeightings.DIRECT_REFERENCE,
            unit=1
        ))
        self.find_or_create(term).add_context(OntologyContextDirectReference())

    def searched_term(self, term: str):
        self.find_or_create(term).add_context(OntologyContextSearched())

    def __iter__(self):
        values = sorted(list(self.term_map.values()), reverse=True)
        return iter(values)

    def as_json(self) -> Any:
        return [value.as_json() for value in self]

    @staticmethod
    def from_search(search_text: str, gene_symbol: Optional[str], selected: Optional[List[str]] = None, server_search: bool = True) -> 'OntologyMatching':
        ontology_matches = OntologyMatching(search_term=search_text, gene_symbol=gene_symbol)
        ontology_matches.populate_relationships()

        if selected:
            for select in selected:
                ontology_matches.select_term(select)

        ONTOLOGY_PATTERN = re.compile("(MONDO|OMIM|HP):[0-9]+")
        if ONTOLOGY_PATTERN.match(search_text):
            ontology_matches.find_or_create(search_text)

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
                        # TODO we would set the label based on the server search previous

                        #label = result.get('label')
                        #if label:
                        #    label = label[0]
                        #
                        #onto = ontology_matches.find_or_create(term_id=o_id)
                        # if not onto.name:
                        #    onto.name = label
                        #onto.add_context(OntologyContextSearched())
                        ontology_matches.searched_term(o_id)
                except:
                    pass
                    # TODO communicate to the user couldn't search mondo text search

        ontology_matches.apply_word_match_score()
        return ontology_matches
