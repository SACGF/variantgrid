import re
from abc import ABC
from dataclasses import dataclass
from re import Match, IGNORECASE
from typing import Optional, TypeVar, Generic, Union, List, Iterable, Type

import django
from django.contrib.auth.models import User
from django.utils.safestring import SafeString

from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild

search_signal = django.dispatch.Signal()
HAS_ALPHA_PATTERN = r"[a-zA-Z]"


@dataclass(frozen=True)
class SearchInput:
    user: User
    search_string: str
    genome_build_preferred: GenomeBuild

    def matches_pattern(self, pattern: str) -> Match:
        return re.match(pattern, self.search_string, IGNORECASE)

    def matches_has_alpha(self) -> bool:
        return bool(self.matches_pattern(HAS_ALPHA_PATTERN))

    def search(self) -> List['SearchResponse']:
        return [response for _, response in search_signal.send(sender=SearchInput, search_input=self) if response.valid_search]


T = TypeVar("T")


class SearchResponseRecordAbstract(ABC, Generic[T]):

    def __init__(self, record: T):
        self.record = record

    @classmethod
    def search_type(cls) -> str:
        pass

    @classmethod
    def from_iterable(cls, iterable: Iterable[T]) -> List:
        return [cls(r) for r in iterable]

    @property
    def genome_build(self) -> Optional[GenomeBuild]:
        return None

    @property
    def annotation_consortia(self) -> Optional[AnnotationConsortium]:
        return None

    # TODO, we might want to add context in so we can get better rendering
    def display(self) -> Union[str, SafeString]:
        return str(self.record)

    def __str__(self):
        return self.display()

    def get_absolute_url(self) -> str:
        return self.record.get_absolute_url()

    @property
    def messages(self) -> Optional[List[str]]:
        return None


class SearchResponse(Generic[T]):

    def __init__(self, response_type: Type[SearchResponseRecordAbstract]):
        self.response_type = response_type
        self.results: List[SearchResponseRecordAbstract] = list()
        self.valid_search = False

    @property
    def search_type(self):
        return self.response_type.search_type()

    def add(self, record: Union[T, SearchResponseRecordAbstract]):
        if not isinstance(record, SearchResponseRecordAbstract):
            record = self.response_type(record)
        self.results.append(self.response_type(record=record))

    def extend(self, iterable: Iterable[T]):
        for r in iterable:
            self.add(r)
        self.mark_valid_search()

    def mark_valid_search(self):
        self.valid_search = True
