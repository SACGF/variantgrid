import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from re import Match, IGNORECASE
from typing import Optional, TypeVar, Generic, Union, List, Iterable, Type, Pattern, Any

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
                    valid_responses.append(response)
                else:
                    # TODO see if there's a more useful way we can pass exceptions?
                    print(caller)
                    report_message("Error during search", 'error', extra_data={"target": str(response), "caller": str(caller)})

        return valid_responses


T = TypeVar("T")


@dataclass
class SearchResult2:
    preview: PreviewData
    genome_build: Optional[GenomeBuild] = None
    annotation_consortium: Optional[AnnotationConsortium] = None
    messages: Optional[List[str]] = None


class SearchResponse:

    def __init__(self, *args):
        self.searched_categories = set()
        self.results = list()
        for arg in args:
            self.add_search_category(arg)

    def add_search_category(self, obj):
        if not isinstance(obj, str):
            obj = obj.preview_category()
        self.searched_categories.add(obj)

    def add_search_result(self, search_result: SearchResult2):
        self.results.append(search_result)
        self.searched_categories.add(search_result.preview.category)

    def add(self, obj: Any, messages: Optional[List[str]] = None, genome_build: Optional[GenomeBuild] = None, annotation_consortium: Optional[GenomeBuild] = None):
        if not isinstance(obj, PreviewData):
            preview = obj.preview
            if not preview:
                raise ValueError(f"{obj} had None preview")
            obj = preview
        if not isinstance(obj, PreviewData):
            raise ValueError(f"Can't add {obj} as search result preview")
        self.add_search_result(SearchResult2(preview=obj, messages=messages, genome_build=genome_build, annotation_consortium=annotation_consortium))

    def extend(self, iterable: Iterable, messages: Optional[List[str]] = None, genome_build: Optional[GenomeBuild] = None, annotation_consortium: Optional[GenomeBuild] = None):
        for obj in iterable:
            self.add(obj, messages=messages, genome_build=genome_build, annotation_consortium=annotation_consortium)
