"""
ConditionResolved: a classification's (or a group's) condition as matched ontology terms and plain text
references, with how multiple terms combine (MultiCondition). Built from the cached
Classification.condition_resolution JSON (ConditionResolvedDict) through ConditionResolved.from_dict and
written back with ConditionResolved.to_json; condition_text_matching.py:apply_condition_resolution is what fills it.
"""
from __future__ import annotations

import operator
from dataclasses import dataclass
from functools import cached_property, reduce
from typing import Optional, TypedDict

import pydantic
from typing_extensions import deprecated

from classification.enums.condition_enums import MultiCondition
from ontology.models import OntologySnake, OntologyTerm, OntologyTermRelation


class ConditionResolvedReferenceDict(TypedDict, total=False):
    term_id: Optional[str]
    name: str  # name of the term if term_id is provided, otherwise from a plain text condition
    count: Optional[int]  # assumed to be 1 if not provided


@deprecated("Use ConditionResolvedReferenceDict")
class ConditionResolvedTermDict(TypedDict, total=False):
    term_id: str
    name: str


class ConditionResolvedDict(TypedDict, total=False):
    """
    Structure of data used to cached resolved condition text again a classification.
    Used to use "resolved_terms" and "plain_text_terms" to store matched OntologyIDs and plain text conditions
    respectively, but now both those are bundled into references.
    Fields are still defined for reading pre-existing JSON.
    """
    display_text: str  # plain text to show to users if not in a position to render links
    sort_text: str  # lower case representation of description
    resolved_terms: list[ConditionResolvedTermDict]
    """
    DEPRECATED, please use references instead
    """

    plain_text_terms: list[str]  # A full list of unresolved plain text conditions
    """
    DEPRECATED, please use references instead
    """

    resolved_join: str

    references: list[ConditionResolvedReferenceDict]  # new method of saving counts


@pydantic.dataclasses.dataclass(frozen=True, config=pydantic.ConfigDict(arbitrary_types_allowed=True))
class ConditionReference:
    """
    Represents an OntologyTerm (just provide the OntologyTerm) or
    Free condition text (just provide name)
    as well as a count of how many times they've appeared.
    Count will be more meaningful when looking at a group of records rather than individual.
    """

    term: Optional[OntologyTerm] = None
    name: Optional[str] = None
    count: int = 1

    @property
    def full_text(self) -> str:
        if term := self.term:
            return f"{term.pk} {self.name}"
        else:
            return self.name

    def to_json(self) -> ConditionResolvedReferenceDict:
        result: ConditionResolvedReferenceDict = {}
        if term := self.term:
            result["term_id"] = term.id
            result["name"] = term.name
        else:
            result["name"] = self.name
        if count := self.count:
            result["count"] = count
        return result

    def __lt__(self, other) -> bool:
        if self.term and other.term:
            return self.term < other.term
        if not self.term and not other.term:
            return self.name < other.name
        if self.term:
            return False
        return True

# TODO: these classes should be reworked as dataclass json

@dataclass(frozen=True)
class ConditionResolved:
    # terms: list[OntologyTerm]
    # plain_text_terms: list[str] = None
    references: list[ConditionReference]
    join: Optional[MultiCondition] = None
    plain_text: Optional[str] = None  # fallback, not populated in all contexts

    def __hash__(self) -> int:
        hash_value = 34543
        if join := self.join:
            hash_value += hash(join)
        if references := self.references:
            hash_value += reduce(operator.add, map(hash, references))
        return hash_value

    @property
    def terms(self) -> list[OntologyTerm]:
        return [t.term for t in self.references if t.term is not None]

    @property
    def plain_text_terms(self) -> list[str]:
        # [] rather than None: callers put this straight into JSON
        return [t.name for t in self.references if t.term is None]

    @staticmethod
    def from_uncounted_terms(
            terms: list[OntologyTerm] = None,
            plain_text_terms: list[str] = None,
            join: Optional[MultiCondition] = None,
            plain_text: Optional[str] = None
    ) -> ConditionResolved:
        references: list[ConditionReference] = []
        if terms:
            references += [ConditionReference(term=term) for term in terms]
        if plain_text_terms:
            references += [ConditionReference(name=plain_text_term) for plain_text_term in plain_text_terms]
        return ConditionResolved(references=references, join=join, plain_text=plain_text)

    @property
    def summary(self) -> str:
        text = ", ".join([term.id for term in self.terms])
        if join := self.join:
            try:
                join = join.value
            except Exception:
                pass
            text = f"{text} {join}"
        return text

    @staticmethod
    def from_dict(condition_dict: ConditionResolvedDict) -> ConditionResolved:
        if not condition_dict:
            return ConditionResolved(references=[])

        join: Optional[MultiCondition] = None
        plain_text = condition_dict.get("display_text")

        if resolved_join_text := condition_dict.get("resolved_join"):
            join = MultiCondition(resolved_join_text)

        if references_list := condition_dict.get("references"):
            references: list[ConditionReference] = []
            reference_dict: ConditionResolvedReferenceDict
            for reference_dict in references_list:
                term: Optional[OntologyTerm] = None
                if term_id := reference_dict.get("term_id"):
                    term = OntologyTerm.get_or_stub_cached(term_id)

                references.append(ConditionReference(
                    term=term,
                    name=reference_dict.get("name"),
                    count=reference_dict.get("count")
                ))
            references.sort()
            return ConditionResolved(references=references, join=join, plain_text=plain_text)
        else:
            # old format with resolved terms and plain_text_terms
            if resolved_terms := condition_dict.get("resolved_terms"):
                terms = list(sorted(OntologyTerm.get_or_stub_cached(term_dict.get("term_id")) for term_dict in resolved_terms))
            else:
                terms = []
            plain_text_terms = condition_dict.get("plain_text_terms")
            return ConditionResolved.from_uncounted_terms(terms=terms, plain_text_terms=plain_text_terms, join=join, plain_text=plain_text)

    @property
    def is_multi_condition(self) -> bool:
        return len(self.references) > 1

    @property
    def single_term(self) -> Optional[OntologyTerm]:
        """
        Only call if not is_multi_condition
        """
        return self.terms[0] if len(self.terms) == 1 else None

    @cached_property
    def mondo_term(self) -> Optional[OntologyTerm]:
        if term := self.single_term:
            return OntologyTermRelation.as_mondo(term)
        else:
            return None

    def as_mondo_if_possible(self) -> ConditionResolved:
        if mondo_term := self.mondo_term:
            return ConditionResolved.from_uncounted_terms(terms=[mondo_term])
        else:
            return self

    def is_same_or_more_specific(self, other: ConditionResolved) -> bool:
        """
        Returns the number of steps to go from this condition to the specific other condition
        Returns None if self doesn't appear to be a descendant of other
        :param other:
        :return:
        """
        if self.is_multi_condition or other.is_multi_condition:
            # when looking at multiple conditions, do not attempt merging unless we're the exact same
            return self.terms == other.terms and self.join == other.join
        elif self.single_term == other.single_term:
            return True
        else:
            if other_mondo := other.mondo_term:
                if self_mondo := self.mondo_term:
                    if other_mondo.index == 1:
                        return 99  # MOND:000001 is always going to be an ancestor
                    descendant_relationships = OntologySnake.check_if_ancestor(descendant=self_mondo,
                                                                               ancestor=other_mondo)
                    return bool(descendant_relationships)

            # terms cant be converted to MONDO and not exact match, just return False
            return False

    @staticmethod
    def more_general_term_if_related(resolved_1: ConditionResolved, resolved_2: ConditionResolved) -> Optional[
        ConditionResolved]:
        more_general: Optional[ConditionResolved] = None
        if resolved_1.is_same_or_more_specific(resolved_2):
            more_general = resolved_2
        elif resolved_2.is_same_or_more_specific(resolved_1):
            more_general = resolved_1

        if more_general:
            # if presented with different types, and we can switch over to MONDO, do so
            if not more_general.is_multi_condition and \
                    resolved_1.single_term.ontology_service != resolved_2.single_term.ontology_service:
                if mondo_term := more_general.mondo_term:
                    more_general = ConditionResolved.from_uncounted_terms(terms=[mondo_term])
            return more_general

        return None

    @staticmethod
    def term_to_dict(term: OntologyTerm) -> ConditionResolvedTermDict:
        term_dict: ConditionResolvedTermDict = {
            "term_id": term.id,
            "name": term.name
        }
        return term_dict

    @property
    def as_plain_text(self) -> str:
        if text := self.plain_text:
            return text
        else:
            return self.to_json()['display_text']

    def to_json(self, include_join: bool = True) -> ConditionResolvedDict:
        jsoned: ConditionResolvedDict
        if self.terms:
            def format_term(term: OntologyTerm) -> str:
                if name := term.name:
                    return f"{term.id} {name}"
                return term.id

            terms = self.terms
            text = ", ".join([format_term(term) for term in terms])
            if self.plain_text_terms:
                text += ",".join(self.plain_text_terms)

            sort_text = ", ".join([term.name for term in terms]).lower()
            join: Optional[MultiCondition] = None
            if len(terms) > 1 and include_join:
                join = self.join or MultiCondition.NOT_DECIDED
                text = f"{text}; {join.label}"

            resolved_term_dicts: list[ConditionResolvedTermDict] = [ConditionResolved.term_to_dict(term) for term in
                                                                    self.terms]
            jsoned: ConditionResolvedDict = {
                "resolved_terms": resolved_term_dicts,
                "resolved_join": join,
                "plain_text_terms": self.plain_text_terms,
                "display_text": text,
                "sort_text": sort_text
            }
            return jsoned
        else:
            jsoned: ConditionResolvedDict = {
                "plain_text_terms": self.plain_text_terms,
                "display_text": ", ".join(pt.lower() for pt in self.plain_text) if self.plain_text else None,
                "sort_text": ", ".join(pt.lower() for pt in self.plain_text) if self.plain_text else None
            }
        return jsoned


    @cached_property
    def join_text(self) -> Optional[str]:
        if join := self.join:
            try:
                return MultiCondition(join).label
            except Exception:
                pass
        return None

    def __eq__(self, other: ConditionResolved) -> bool:
        if (s_terms := self.terms) and (o_terms := other.terms):
            return s_terms == o_terms and self.join == other.join
        elif self.terms or other.terms:
            return False
        else:
            return self.plain_text == other.plain_text

    def __lt__(self, other: ConditionResolved) -> bool:
        self_terms = self.terms or []
        other_terms = other.terms or []
        if len(self_terms) != len(other_terms):
            return len(self_terms) < len(other_terms)
        if len(self_terms) >= 1:
            return self_terms[0] < other_terms[0]
        return (self.plain_text or '') < (other.plain_text or '')
