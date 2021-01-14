import operator
import re
from dataclasses import dataclass
from functools import reduce
from typing import Dict, Optional, List, Any, Set

import requests
from django.db.models import Q
from django.db.models.functions import Length

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


class OntologyContextSearched(OntologyContext):
    """
    Indicates that this was found from doing a text search
    e.g. search for "rett syndrome" and found "MONDO:0010726" based on text match.
    """

    @property
    def name(self):
        return "searched"

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
        return "referenced"

    def as_json(self) -> Any:
        return True


class OntologyContextSelected(OntologyContext):
    """
    Indicates that this value has already been selected for the text
    """

    @property
    def name(self):
        return "selected"

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
        return "Monarch"

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
        return "PanelApp AU"

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
        return "parent_of_gene_relation"

    def as_json(self):
        return {
            "matching_children": self.children,
        }


@dataclass
class OntologyMatchScorePart:
    description: str
    score: float

    def as_json(self):
        return {
            "description": self.description,
            "score": self.score
        }


class OntologyMeta:

    def __init__(self, term_id: str):
        self.term = OntologyTerm.get_or_stub(term_id)
        self.parent: Optional[OntologyTerm] = None
        self.scores: List[OntologyMatchScorePart] = list()

        if not self.term.is_stub:
            # FIXME multiple parents are valid
            self.parent = OntologyTermRelation.parents_of(self.term).first()
        self.contexts: Dict[str, OntologyContext] = dict()

    @property
    def score(self) -> float:
        score = 0
        for score_part in self.scores:
            score += score_part.score
        return score

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
            "url": self.term.url,
            "title": self.term.name,
            "definition": self.term.definition,
            "contexts": contexts,
            "scores": scores,
            "score": self.score
        }

        if self.parent:
            data["parent_id"] = self.parent.id
            data["parent_title"] = self.parent.name
            """
            sibling_count = MonarchDiseaseOntologyRelationship.objects.filter(object=self.parent, relationship="is_a").count() - 1
            child_count = MonarchDiseaseOntologyRelationship.objects.filter(object=self.object, relationship="is_a").count()

            data["sibling_count"] = sibling_count
            data["children_count"] = child_count
            """
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
                mondo.add_score(OntologyMatchScorePart("Base score within gene context", 50.0))
            else:
                mondo.add_score(OntologyMatchScorePart("Base score outside of gene context", 100.0))

        return mondo

    def apply_word_match_score(self):
        for onto in self.term_map.values():
            self.add_word_match_score(onto)

    def add_word_match_score(self, mondo: OntologyMeta):
        if search_terms := self.search_terms:
            if mondo.term.is_stub:
                mondo.add_score(OntologyMatchScorePart(
                    "No local copy of this term",
                    -100
                ))
            else:
                match_terms = set(tokenize_condition_text(normalize_condition_text(mondo.term.name))) - SKIP_TERMS
                superfluous_words = match_terms.difference(search_terms)
                missing_words = search_terms.difference(match_terms)

                missing_word_ratio = 0 if not match_terms else (float(len(missing_words)) / float(len(match_terms)))
                superfluous_word_ratio =  0 if not match_terms else (float(len(superfluous_words)) / float(len(match_terms)))

                if missing_word_ratio:
                    mondo.add_score(OntologyMatchScorePart(
                        f"Missing words {pretty_set(missing_words)} @ ratio {missing_word_ratio:.2f} * 100",
                        missing_word_ratio * -100.0
                    ))

                if superfluous_word_ratio:
                    mondo.add_score(OntologyMatchScorePart(
                        f"Superfluous words {pretty_set(superfluous_words)} @ ratio {superfluous_word_ratio:.2f} * 50",
                        superfluous_word_ratio * -50.0
                    ))

                if superfluous_words and not self.gene_symbol:
                    for word in superfluous_words:
                        if SUB_TYPE.match(word):
                            mondo.add_score(OntologyMatchScorePart(
                                "Superfluous words includes subtype \"{word}\"",
                                20.0
                            ))

    def populate_relationships(self):
        if gene_symbol := self.gene_symbol:
            update_gene_relations(gene_symbol=gene_symbol)
            snakes = OntologySnake.terms_for_gene_symbol(gene_symbol=gene_symbol, desired_ontology=OntologyService.MONDO)
            for snake in snakes:
                # always convert to MONDO for now
                mondo_term = snake.leaf_term
                mondo_meta = self.find_or_create(mondo_term.id)
                gene_relation = snake.paths[0]
                gene_relationshpis_via = list()
                if gene_relation.relation == OntologyRelation.PANEL_APP_AU:
                    mondo_meta.add_context(OntologyContextPanelApp(
                        gene_symbol=gene_symbol,
                        omim_id=0, # TODO populate the OMIM ID
                        phenotype_row=gene_relation.extra.get("phenotype_row"),
                        evidence=gene_relation.extra.get("evidence")
                    ))
                    gene_relationshpis_via.append("PanelApp AU")
                else:
                    mondo_meta.add_context(OntologyContextMonarchLink(
                        gene_symbol=gene_symbol,
                        relation=gene_relation.relation
                    ))
                    gene_relationshpis_via.append("Other")

                gene_relationshpis_via_str = ", ".join(gene_relationshpis_via)
                mondo_meta.add_score(OntologyMatchScorePart(
                    f"Established gene relationship (through {gene_relationshpis_via_str})", 50.0
                ))

                for parent_term in OntologyTermRelation.parents_of(mondo_term):
                    parent_meta = self.find_or_create(parent_term.id)
                    parent_meta.add_context(OntologContextChildGeneRelationship(
                        child_term_id=mondo_term.id
                    ))
                    # FIXME don't rely on magic string for the name of the contexts
                    if "Monarch" not in mondo_meta.contexts and "PanelApp AU" not in mondo_meta.contexts:
                        mondo_meta.add_score(OntologyMatchScorePart(
                            "Parent of term with established gene relationship", 25.0
                        ))

    def select_term(self, term: str):
        self.find_or_create(term).add_context(OntologyContextSelected())

    def reference_term(self, term: str):
        onto_meta = self.find_or_create(term)
        onto_meta.add_score(OntologyMatchScorePart(
            "Directly referenced", 1000.0
        ))
        self.find_or_create(term).add_context(OntologyContextDirectReference())

    def searched_term(self, term: str):
        self.find_or_create(term).add_context(OntologyContextSearched())

    def as_json(self) -> Any:
        return [value.as_json() for value in self.term_map.values()]

    def __iter__(self):
        # TODO probably should sort as well
        return iter(self.term_map.values())

    @staticmethod
    def from_search(search_text: str, gene_symbol: Optional[str], selected: Optional[List[str]] = None, server_search: bool = True):
        ontology_matches = OntologyMatching(search_term=search_text, gene_symbol=gene_symbol)
        ontology_matches.populate_relationships()

        if selected:
            for select in selected:
                ontology_matches.select_term(select)

        # TODO change to any pattern, not just MONDO
        MONDO_PATTERN = re.compile("MONDO:[0-9]+")
        if MONDO_PATTERN.match(search_text):
            ontology_matches.find_or_create(search_text)

        else:
            # Client search
            # this currently requires all words to be present
            search_terms = set(tokenize_condition_text(search_text))
            qs = OntologyTerm.objects.filter(ontology_service=OntologyService.MONDO)
            qs = qs.filter(reduce(operator.and_, [Q(name__icontains=term) for term in search_terms]))
            qs = qs.order_by(Length('name')).values_list("id", flat=True)
            for result in qs[0:20]:
                ontology_matches.searched_term(result)

            # ALSO, do client side search
            if server_search:
                try:
                    results = requests.get(f'https://api.monarchinitiative.org/api/search/entity/autocomplete/{search_text}', {
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
