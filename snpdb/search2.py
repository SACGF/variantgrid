import re
from dataclasses import dataclass, field
from re import Match, IGNORECASE
from typing import Optional, TypeVar, Union, List, Iterable, Pattern, Any, Set
import django
from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet
from genes.models_enums import AnnotationConsortium
from library.log_utils import report_message
from library.preview_request import PreviewData
from snpdb.models import GenomeBuild, Variant

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

    @property
    def search_words(self) -> List[str]:
        words = [word.strip() for word in self.search_string.split(" ")]
        return [word for word in words if word]

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

    def get_visible_variants(self, genome_build: GenomeBuild) -> QuerySet[Variant]:
        """ Shariant wants to restrict search to only classified variants """

        from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
        from annotation.models import AnnotationVersion
        from classification.models import Classification

        annotation_version = AnnotationVersion.latest(genome_build)
        variant_qs = get_variant_queryset_for_annotation_version(annotation_version)
        variant_qs = variant_qs.filter(Variant.get_contigs_q(genome_build))  # restrict to build
        if settings.SEARCH_VARIANT_REQUIRE_CLASSIFICATION_FOR_NON_ADMIN and not self.user.is_superuser:
            variant_qs = variant_qs.filter(Classification.get_variant_q(self.user, genome_build))

        return variant_qs

    @property
    def genome_builds(self) -> Set[GenomeBuild]:
        return GenomeBuild.builds_with_annotation()


T = TypeVar("T")


@dataclass
class SearchResult2:
    preview: PreviewData
    genome_builds: Optional[Set[GenomeBuild]] = None
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

    def add(self, obj: Any, messages: Optional[List[str]] = None, genome_builds: Optional[Set[GenomeBuild]] = None, annotation_consortium: Optional[GenomeBuild] = None):
        if not isinstance(obj, PreviewData):
            preview = obj.preview
            if not preview:
                raise ValueError(f"{obj} had None preview")
            obj = preview
        if not isinstance(obj, PreviewData):
            raise ValueError(f"Can't add {obj} as search result preview")
        self.add_search_result(SearchResult2(preview=obj, messages=messages, genome_builds=genome_builds, annotation_consortium=annotation_consortium))

    def extend(self, iterable: Iterable, messages: Optional[List[str]] = None, genome_builds: Optional[Set[GenomeBuild]] = None, annotation_consortium: Optional[GenomeBuild] = None):
        for obj in iterable:
            self.add(obj, messages=messages, genome_builds=genome_builds, annotation_consortium=annotation_consortium)
