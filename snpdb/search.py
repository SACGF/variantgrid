import operator
import re
import itertools
from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from functools import cached_property, reduce
from re import IGNORECASE
from typing import List, Set, Optional, Type, Pattern, Callable, Any, Match, Union, Dict, Tuple

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet, Q
from django.dispatch import Signal

from library.enums.log_level import LogLevel
from library.log_utils import report_exc_info, report_message
from library.preview_request import PreviewModelMixin, PreviewData
from library.utils import clean_string
from snpdb.models import UserSettings, GenomeBuild, Variant, Allele


search_signal = Signal()
HAS_ALPHA_PATTERN = re.compile(r"[a-zA-Z]")
HAS_ANYTHING = re.compile(r".")
HAS_3_ALPHA_MIN = re.compile(r"[a-zA-Z]{3,}")
HAS_3_ANY = re.compile(r"\S{3,}")
_SPLIT_GAPS = re.compile(r"[\s,]+")


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
        return [word for word in _SPLIT_GAPS.split(self.search_string) if word]

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

    @cached_property
    def genome_builds(self) -> Set[GenomeBuild]:
        return set(GenomeBuild.builds_with_annotation().all())


@dataclass(frozen=True)
class SearchInputInstance:
    expected_type: Type
    search_input: SearchInput
    match: Match

    @property
    def classify(self) -> bool:
        return self.search_input.classify

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
class SearchMessage:
    message: str
    severity: LogLevel = LogLevel.WARNING


@dataclass
class SearchMessagesForType:
    search_message: SearchMessage
    search_info: Any

    @property
    def severity(self):
        return self.search_message.severity

    @property
    def message(self):
        return self.search_message.message


class SearchResultMatchStrength(int, Enum):
    FUZZY_MATCH = 1
    STRONG_MATCH = 2
    ID_MATCH = 3


@dataclass
class SearchResult:
    preview: PreviewData
    messages: Optional[List[str]] = None
    sub_name: Optional[str] = None
    match_strength: SearchResultMatchStrength = None
    ignore_genome_build_mismatch: bool = False

    @property
    def effective_match_strength(self):
        return self.match_strength or SearchResultMatchStrength.STRONG_MATCH

    def _preview_icon_severity(self, severity: str):
        icon = self.preview.icon
        if not severity:
            return icon

        # icon has a background colour built in, can replace that
        if 'bg-' in icon:
            return re.sub(r'bg-\w+', f'bg-{severity}', icon)
        else:
            return f"{icon} text-{severity}"

    @property
    def preview_icon_with_severity(self):
        if self.messages:
            return self._preview_icon_severity('warning')
        else:
            return self._preview_icon_severity('success')

    @property
    def search_type(self):
        return self.preview.category

    @property
    def genome_builds(self):
        return self.preview.genome_builds

    @property
    def genome_build_names(self):
        if gbs := self.preview.genome_builds:
            return list(sorted(gbs))
        return []

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


def _default_make_search_result(obj) -> Optional['SearchResult']:
    if isinstance(obj, SearchResult):
        return obj
    else:
        return SearchResult(preview=obj.preview)


@dataclass
class _SearchResultFactory:
    source: Any
    converter: Callable[[Any], SearchResult]
    extra_messages: List[str]

    @staticmethod
    def convert(output: Any):
        source: Any
        factory: Callable[[Any], SearchResult] = _default_make_search_result
        extra_messages: List[str] = []
        if isinstance(output, (tuple, list)):
            source = output[0]
            for extra in output[1:]:
                if isinstance(extra, str):
                    extra_messages.append(extra)
                elif isinstance(extra, list):
                    extra_messages += extra
                elif isinstance(extra, Callable):
                    factory = extra
                else:
                    raise ValueError(f"Search method yielded extra {extra}, expected str, list[str] or callable")
        else:
            source = output

        return _SearchResultFactory(source=source, converter=factory, extra_messages=extra_messages)

    @property
    def _iterable(self):
        if isinstance(self.source, (list, QuerySet)):
            return self.source
        else:
            return [self.source]

    def __len__(self) -> int:
        if isinstance(self.source, list):
            return len(self.source)
        elif isinstance(self.source, QuerySet):
            return self.source.count()
        else:
            return 1

    def iterate(self, limit: int = 20):
        for obj in self._iterable[:limit]:
            search_result = self.converter(obj)
            if search_result is None:
                continue
            elif not isinstance(search_result, SearchResult):
                raise ValueError(f"search result factory returned {search_result} instead of a SearchResult")
            if self.extra_messages:
                search_result.messages += self.extra_messages
            yield search_result


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
    results: List[SearchResult] = field(default_factory=list)
    messages: List[SearchMessage] = field(default_factory=list)
    total_count: int = 0

    @property
    def search_messages_with_type(self):
        return [SearchMessagesForType(search_message=message, search_info=self) for message in self.messages]

    @property
    def preview_category(self):
        # django templates don't handle {{ UserPreview.preview_category }} for some unknown reason, works for all models classes
        return self.search_type_effective.preview_category()

    @property
    def preview_icon(self):
        return self.search_type_effective.preview_icon()

    def add_search_result(self, search_result: SearchResult):
        self.results.append(search_result)

    @property
    def search_type_effective(self) -> PreviewModelMixin:
        if settings.PREFER_ALLELE_LINKS and self.search_type == Variant:
            return Allele
        return self.search_type

    @property
    def optimized_results(self):
        # used to just test if self.search_type == Variant, but DJango in debug mode seems to get confused by that on dev restarts
        if self.results and settings.PREFER_ALLELE_LINKS and self.search_type.preview_category() == "Variant":
            allele_to_variants: Dict[Allele, List[SearchResult]] = defaultdict(list)
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
                        SearchResult(preview=allele_preview, messages=messages)
                    )
            return regular_results
        return self.results


@dataclass
class SearchCount:
    category: str
    resolved: int = 0
    total: int = 0

    @property
    def excluded_results(self) -> bool:
        return self.total > self.resolved


class CombinedSearchResponses:

    def __init__(self, search_input: SearchInput, responses: List[SearchResponse]):
        self.search_input = search_input
        self.responses = responses

    def search_counts(self) -> List[SearchCount]:
        counts: Dict[str, SearchCount] = {}
        for response in self.responses:
            # FIXME if the responses expected return type doesn't match a SearchResult's actual category
            # these counts all go to bunk, probably best to disallow that
            category = response.search_type.preview_category()
            sc = counts.get(category)
            if sc is None:
                sc = SearchCount(category=category)
                counts[category] = sc
            sc.resolved += len(response.results)
            sc.total += response.total_count

        return list(sorted((sc for sc in counts.values() if sc.resolved), key=lambda x: x.category))

    @cached_property
    def search_types(self) -> Set[str]:
        search_types_set = set()
        for response in self.responses:
            if response.matched_pattern:
                search_types_set.add(response.search_type.preview_category())
        return search_types_set

    @cached_property
    def ordered_responses(self):
        return list(sorted(self.responses, key=lambda x: x.search_type_effective.preview_category()))

    @cached_property
    def results(self) -> List[SearchResult]:
        result_list = list(itertools.chain.from_iterable([response.optimized_results for response in self.ordered_responses]))
        # TODO, sort these
        return result_list

    @cached_property
    def messages(self):
        message_list = list(itertools.chain.from_iterable([response.search_messages_with_type for response in self.responses]))
        # TODO, sort these
        return message_list

    def single_preferred_result(self):
        if not self.messages and len(self.results) == 1:
            first_result = self.results[0]
            if not first_result.messages:
                return first_result

    @property
    def summary(self) -> str:
        return f"{self.search_input.search_string}' calculated {len(self.results)} results."

    @cached_property
    def grouped_search_infos(self) -> List[SearchResponse]:
        return sorted(self.responses, key=lambda sr: (sr.preview_category, sr.sub_name.upper() if sr.sub_name else ""))


def search_data(user: User, search_string: str, classify: bool = False) -> CombinedSearchResponses:
    user_settings = UserSettings.get_for_user(user)
    search_input = SearchInput(
        user=user,
        search_string=clean_string(search_string),
        genome_build_preferred=user_settings.default_genome_build,
        classify=classify
    )
    responses = search_input.search()
    return CombinedSearchResponses(search_input, responses)


def search_receiver(
        search_type: Optional[Type[PreviewModelMixin]],
        pattern: Pattern = HAS_ANYTHING,
        admin_only: bool = False,
        sub_name: Optional[str] = None,
        example: Optional[SearchExample] = None,
        match_strength: Optional[SearchResultMatchStrength] = SearchResultMatchStrength.STRONG_MATCH
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
    :param match_strength: Provide an overall match strength to be used if the individual search results don't set it
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
                    # as Variants get merged into Alleles, we want to avoid limiting them (except under extreme conditions)
                    limit = 1000 if search_type.preview_category() == "Variant" else 50
                    total_count = 0
                    for result in func(SearchInputInstance(expected_type=search_type, search_input=search_input, match=match)):
                        if result is None:
                            raise ValueError(f"Search {sender.__name__} returned None")
                        elif isinstance(result, SearchMessage):
                            response.messages.append(result)
                        else:

                            # need to make sure Variable is allowed to return more results as they get halved into alleles
                            factory = _SearchResultFactory.convert(result)

                            total_count += len(factory)

                            if limit > 0:
                                for search_result in factory.iterate(limit=limit):
                                    if match_strength and not search_result.match_strength:
                                        search_result.match_strength = match_strength
                                    if sub_name:
                                        search_result.sub_name = sub_name

                                    response.add_search_result(search_result)
                                    limit -= 1
                    response.total_count = total_count

            except Exception as e:
                # TODO, determine if the Exception type is valid for users or not
                response.messages.append(SearchMessage(str(e), severity=LogLevel.ERROR))
                print(f"Error handling search_receiver on {func}")
                report_exc_info()

            return response

        search_signal.connect(search_func)
        return search_func

    return _decorator
