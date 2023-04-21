import operator
import re
from abc import ABC
from collections import defaultdict
from dataclasses import dataclass, field
from functools import reduce, cached_property
from re import Match, IGNORECASE
from typing import Optional, TypeVar, Union, List, Iterable, Pattern, Any, Set, Callable, Type, Dict
from django.dispatch import Signal
from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet, Q
from genes.models_enums import AnnotationConsortium
from library.log_utils import report_message, report_exc_info
from library.preview_request import PreviewData, PreviewModelMixin
from snpdb.models import GenomeBuild, Variant, Allele
from variantgrid.perm_path import get_visible_url_names

search_signal = Signal()
HAS_ALPHA_PATTERN = re.compile(r"[a-zA-Z]")
HAS_ANYTHING = re.compile(r".")


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
class SearchResult2:
    preview: PreviewData
    messages: Optional[List[str]] = None

    @property
    def genome_builds(self):
        if gb := self.preview.genome_build:
            return {gb}
        return None

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


@dataclass(frozen=True)
class SearchExample:
    note: Optional[str] = None
    example: Optional[str] = None


@dataclass
class SearchResponse:
    search_input: SearchInput
    search_type: PreviewModelMixin
    matched_pattern: bool = False
    sub_name: Optional[str] = None
    admin_only: bool = False
    example: Optional[SearchExample] = None
    results: List[SearchResult2] = field(default_factory=list)

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
        if settings.PREFER_ALLELE_LINKS:
            if self.search_type == Variant:
                allele_to_variants: Dict[Allele, List[SearchResult2]] = defaultdict(list)
                regular_results = list()
                had_allele = False
                for result in self.results:
                    if obj := result.preview.obj:
                        if isinstance(obj, Variant):
                            if allele := obj.allele:
                                allele_to_variants[allele].append(result)
                                had_allele = True
                    if not had_allele:
                        regular_results.append(result)

                if allele_to_variants:
                    for allele, variant_results in allele_to_variants.items():
                        genome_builds = set([vr.preview.genome_build for vr in variant_results])
                        has_preferred_allele = self.search_input.genome_build_preferred in genome_builds

                        # FIXME we already have support for genome build ranking, change to use that
                        messages = reduce(operator.__add__, (variant.messages for variant in variant_results if variant.messages), [])
                        if not has_preferred_allele:
                            genome_build_str = ", ".join([gb.name for gb in sorted(genome_builds)])
                            messages.append(f"Found in {genome_build_str}")

                        regular_results.append(
                            # FIXME remove redundant messages
                            SearchResult2(preview=allele.preview, messages=messages)
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
    :param admin_only: Is this search for admin users only, wont be invokved or displayed for non admins.
    :param sub_name: Are there multiple search implementations for the same search_type, if so give them a sub_name
    :param example: An example (to be presented to the user) of how to use this search
    :return: A wrapped call that will obey the above (e.g. preview_enabled(), admin_only) and return a SearchResponse
    """
    def _decorator(func):
        def search_func(sender: Any, search_input: SearchInput, **kwargs):
            try:
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

                if match := pattern.search(search_input.search_string):
                    response.matched_pattern = True
                    for result in func(SearchInputInstance(expected_type=search_type, search_input=search_input, match=match)):
                        for search_result in SearchResult2.convert(result):
                            response.add_search_result(search_result)

                return response
            except:
                # FIXME make this an overall search error
                print(f"Error handling search_receiver on {func}")
                report_exc_info()
                raise

        search_signal.connect(search_func)
        return search_func

    return _decorator