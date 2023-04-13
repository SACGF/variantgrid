import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from re import Match, IGNORECASE
from typing import Optional, TypeVar, Generic, Union, List, Iterable, Type, Pattern

import django
from django.contrib.auth.models import User
from django.utils.safestring import SafeString

from genes.models_enums import AnnotationConsortium
from library.log_utils import report_message
from library.preview_request import PreviewData
from library.utils import pretty_label
from snpdb.models import GenomeBuild

search_signal = django.dispatch.Signal()
HAS_ALPHA_PATTERN = r"[a-zA-Z]"


@dataclass(frozen=True)
class SearchInput:
    user: User
    search_string: str
    genome_build_preferred: GenomeBuild

    def matches_pattern(self, pattern: Union[str, Pattern]) -> Match:
        if isinstance(pattern, str):
            return re.match(pattern, self.search_string, IGNORECASE)
        else:
            return pattern.match(self.search_string)

    def matches_has_alpha(self) -> bool:
        return bool(self.matches_pattern(HAS_ALPHA_PATTERN))

    def search(self) -> List['SearchResponse']:
        valid_responses: List[SearchResponse] = []
        response_tuples = search_signal.send_robust(sender=SearchInput, search_input=self)
        for caller, response in response_tuples:
            if response:
                if isinstance(response, SearchResponse):
                    if response.valid_search:
                        valid_responses.append(response)
                else:
                    # TODO see if there's a more useful way we can pass exceptions?
                    report_message("Error during search", 'error', extra_data={"target": str(response), "caller": str(caller)})

        return valid_responses


T = TypeVar("T")


class SearchResponseRecordAbstract(ABC, Generic[T]):

    def __init__(self, record: T, search_input: Optional[SearchInput] = None):
        self.record = record
        self.search_input = search_input

    @classmethod
    @abstractmethod
    def result_class(cls) -> Type[T]:
        pass

    @classmethod
    def category(cls) -> str:
        return pretty_label(cls.result_class()._meta.verbose_name)

    @property
    def preview(self) -> PreviewData:
        return self.record.preview

    @classmethod
    def from_iterable(cls, iterable: Iterable[T]) -> List:
        return [cls(r) for r in iterable]

    @property
    def genome_build(self) -> Optional[GenomeBuild]:
        return None

    @property
    def annotation_consortia(self) -> Optional[AnnotationConsortium]:
        return None

    @property
    def messages(self) -> Optional[List[str]]:
        return None


class SearchResponse(Generic[T]):

    def __init__(self, response_type: Type[SearchResponseRecordAbstract], search_input: Optional[SearchInput] = None):
        self.response_type = response_type
        self.results: List[SearchResponseRecordAbstract] = []
        self.valid_search = False
        self.search_input = search_input

    @property
    def search_type(self):
        return self.response_type.category()

    def add(self, record: Union[T, SearchResponseRecordAbstract]):
        if not isinstance(record, SearchResponseRecordAbstract):
            record = self.response_type(record, search_input=self.search_input)
        self.results.append(record)

    def extend(self, iterable: Iterable[T]):
        for r in iterable:
            self.add(r)
        self.mark_valid_search()

    def mark_valid_search(self):
        self.valid_search = True
