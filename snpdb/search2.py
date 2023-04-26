import operator
import re
from collections import defaultdict
from dataclasses import dataclass, field
from functools import reduce, cached_property
from re import Match, IGNORECASE
from typing import Optional, Union, List, Pattern, Any, Set, Callable, Type, Dict
from django.dispatch import Signal
from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet, Q
from library.log_utils import report_message, report_exc_info
from library.preview_request import PreviewData, PreviewModelMixin
from snpdb.models import GenomeBuild, Variant, Allele

search_signal = Signal()
HAS_ALPHA_PATTERN = re.compile(r"[a-zA-Z]")
HAS_ANYTHING = re.compile(r".")


@dataclass(frozen=True)
class SearchInput:
    user: User
    search_string: str
    genome_build_preferred: GenomeBuild
    classify: bool = False

    def matches_pattern(self, pattern: Union[str, Pattern]) -> Match:
        if isinstance(pattern, str):
            return re.match(pattern, self.search_string, IGNORECASE)
        else:
            return pattern.match(self.search_string)

    @property
    def search_words(self) -> List[str]:
        words = [word.strip() for word in self.search_string.split(" ")]
        return [word for word in words if word]

    def q_words(self, field_name: str = "name", test: str = "icontains") -> Q:
        """
        Returns Qs and-ed together to split search into words and make sure the words individually are contained
        in the result. So if you searched for "one two" you would get records with "one thing and two" because the
        words are contained, just not necessarily in the order given.
        :param field_name: name by default, The field that should be checked against e.g. "name" results in "name__icontains=<word>"
        :param test: icontains by default, alternatively can do contains.
        :return: A Q to filter your query set with
        """
        if words := self.search_words:
            qs: List[Q] = []
            for word in words:
                qs.append(Q(**{f"{field_name}__{test}": word}))
            return reduce(operator.and_, qs)
        else:
            raise ValueError("No tokens found in search, can't generate q_words")

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


@dataclass(frozen=True)
class SearchInputInstance:
    expected_type: Type
    search_input: SearchInput
    match: Match

    # def matches_pattern(self, pattern: Union[str, Pattern]) -> Match:
    #     return self.search_input.matches_pattern(pattern)

    @cached_property
    def results(self):
        return SearchResponse(self.expected_type)

    @property
    def user(self):
        return self.search_input.user

    @property
    def search_string(self):
        return self.search_input.search_string

    @property
    def genome_builds(self) -> Set[GenomeBuild]:
        return self.search_input.genome_builds

    def get_visible_variants(self, genome_build: GenomeBuild) -> QuerySet[Variant]:
        return self.search_input.get_visible_variants(genome_build=genome_build)

    @property
    def search_words(self) -> List[str]:
        return self.search_input.search_words

    def q_words(self, field_name: str = "name", test: str = "icontains") -> Q:
        return self.search_input.q_words(field_name=field_name, test=test)


@dataclass
class SearchWarning:
    message: str


@dataclass
class SearchWarningForType:
    warning: SearchWarning
    search_info: Any

    @property
    def message(self):
        return self.warning.message


@dataclass
class SearchResult2:
    preview: PreviewData
    messages: Optional[List[str]] = None

    @property
    def search_type(self):
        # Just here to bridge the gap between SeearchResult and SearchResult2
        return self.preview.category

    @property
    def genome_builds(self):
        return self.preview.genome_builds

    @property
    def annotation_consortium(self):
        if ac := self.preview.annotation_consortium:
            return ac
        return None

    def add_message(self, message: str):
        if self.messages is None:
            self.messages = [message]
        else:
            self.messages.append(message)

    @staticmethod
    def convert(data: Any) -> List['SearchResult2']:
        if isinstance(data, SearchResult2):
            return [data]
        if isinstance(data, (tuple, list)):
            data = list(data)
            messages = []
            for entry in data[1:]:
                if isinstance(entry, list):
                    messages.extend(entry)
                elif isinstance(entry, str):
                    messages.append(entry)
                elif entry is None:
                    pass
                else:
                    raise ValueError(f"Not sure what to do with search result extra {entry}")
            previews = SearchResult2._convert_to_previews(data[0])
            return [SearchResult2(preview=preview, messages=messages) for preview in previews]
        else:
            previews = SearchResult2._convert_to_previews(data)
            return [SearchResult2(preview=preview) for preview in previews]

    @staticmethod
    def _convert_to_previews(data: Any) -> List[PreviewData]:
        if isinstance(data, PreviewData):
            return [data]
        elif isinstance(data, QuerySet):
            return list(obj.preview for obj in data)
        else:
            return [data.preview]

    @property
    def header_extra(self) -> List[str]:
        extras = []
        if ac := self.preview.annotation_consortium:
            extras.append(ac.name)
        if self.preview.genome_builds:
            extras.extend((genome_build.name for genome_build in sorted(self.genome_builds)))
        return extras


@dataclass(frozen=True)
class SearchExample:
    note: Optional[str] = None
    examples: Optional[List[str]] = None


@dataclass
class SearchResponse:
    search_input: SearchInput
    search_type: PreviewModelMixin
    matched_pattern: bool = False
    sub_name: Optional[str] = None
    admin_only: bool = False
    example: Optional[SearchExample] = None
    results: List[SearchResult2] = field(default_factory=list)
    warnings: List[SearchWarning] = field(default_factory=list)

    @property
    def search_warnings_with_type(self):
        return [SearchWarningForType(warning=warning, search_info=self) for warning in self.warnings]

    @property
    def preview_category(self):
        # django templates don't handle {{ UserPreview.preview_category }} for some unknown reason, works for all models classes
        return self.search_type_effective.preview_category()

    @property
    def preview_icon(self):
        return self.search_type_effective.preview_icon()

    def add_search_result(self, search_result: SearchResult2):
        self.results.append(search_result)

    @property
    def search_type_effective(self) -> PreviewModelMixin:
        if settings.PREFER_ALLELE_LINKS and self.search_type == Variant:
            return Allele
        return self.search_type

    @property
    def optimized_results(self):
        if self.results and settings.PREFER_ALLELE_LINKS and self.search_type == Variant:
            allele_to_variants: Dict[Allele, List[SearchResult2]] = defaultdict(list)
            regular_results = list()
            for result in self.results:
                had_allele = False
                if obj := result.preview.obj:
                    if isinstance(obj, Variant):
                        if allele := obj.allele:
                            allele_to_variants[allele].append(result)
                            had_allele = True
                if not had_allele:
                    regular_results.append(result)

            if allele_to_variants:
                for allele, variant_results in allele_to_variants.items():
                    genome_builds = set().union(*[vr.preview.genome_builds for vr in variant_results if vr.preview.genome_builds])

                    all_messages = set().union(*(set(variant.messages) for variant in variant_results if variant.messages))
                    messages = sorted(all_messages)

                    allele_preview = allele.preview
                    allele_preview.genome_builds = genome_builds

                    regular_results.append(
                        SearchResult2(preview=allele_preview, messages=messages)
                    )
            return regular_results
        return self.results


def search_receiver(
        search_type: Optional[Type[PreviewModelMixin]] = None,
        pattern: Pattern = HAS_ANYTHING,
        admin_only: bool = False,
        sub_name: Optional[str] = None,
        example: Optional[SearchExample] = None
    ) -> Callable[[Callable], Callable[[SearchInput], SearchResponse]]:
    """
    Wrap around a method that takes a SearchInputInstance and yields QuerySets and/or individual objects, optionally as part of a tuple
    where the subsequent items are validation messages.
    Note that the wrapped method will take a SearchInput and return a SearchResponse object.
    :param search_type: Something that has preview_category(), preview_icon() and preview_enabled() for overall properties of the search
    :param pattern: A regex pattern that must be met for the search to be invoked. The search will be passed a SearchInputInstance with the match result
    :param admin_only: Is this search for admin users only, won't be invoked or displayed for non admins.
    :param sub_name: Are there multiple search implementations for the same search_type, if so give them a sub_name
    :param example: An example (to be presented to the user) of how to use this search
    :return: A wrapped call that will obey the above (e.g. preview_enabled(), admin_only) and return a SearchResponse
    """
    def _decorator(func):
        def search_func(sender: Any, search_input: SearchInput, **kwargs):

            if admin_only and not search_input.user.is_superuser:
                # if only for admin users, don't include it for non admin
                return None
            if search_type and not search_type.preview_enabled():
                #
                return None

            response = SearchResponse(
                search_input=search_input,
                search_type=search_type,
                admin_only=admin_only,
                sub_name=sub_name,
                example=example)

            try:
                if match := pattern.search(search_input.search_string):
                    response.matched_pattern = True
                    for result in func(SearchInputInstance(expected_type=search_type, search_input=search_input, match=match)):
                        if result is None:
                            raise ValueError(f"Search {sender.__name__} returned None")
                        if isinstance(result, SearchWarning):
                            response.warnings.append(result)
                        for search_result in SearchResult2.convert(result):
                            response.add_search_result(search_result)
            except Exception as e:
                response.warnings.append(SearchWarning(str(e)))
                print(f"Error handling search_receiver on {func}")
                report_exc_info()

            return response

        search_signal.connect(search_func)
        return search_func

    return _decorator
