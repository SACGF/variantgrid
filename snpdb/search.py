import operator
import re
import itertools
from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from functools import cached_property, reduce
from re import IGNORECASE
from typing import List, Set, Optional, Type, Pattern, Callable, Any, Match, Union, Dict

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet, Q
from django.dispatch import Signal
from more_itertools import take

from library.enums.log_level import LogLevel
from library.log_utils import report_exc_info, report_message, log_level_to_int, log_level_to_bootstrap
from library.preview_request import PreviewCoordinator, PreviewData
from library.utils import clean_string, first
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
                    # note this doesn't happen as exceptions during search are handled by the search_receiver
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


@dataclass(frozen=True)
class SearchMessageOverall:
    message: str
    severity: LogLevel = LogLevel.WARNING

    @property
    def _sort_order(self):
        return log_level_to_int(self.severity), self.message

    def __lt__(self, other) -> bool:
        return self._sort_order < other._sort_order

    def __post_init__(self):
        if not isinstance(self.message, str):
            raise ValueError(f"Created a Search Message with something other than string {self.message}")


@dataclass
class SearchMessagesOverallForType:
    search_message: SearchMessageOverall
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


@dataclass(frozen=True)
class SearchMessage:
    message: str
    severity: LogLevel = LogLevel.WARNING
    genome_build: Optional[GenomeBuild] = None

    def __post_init__(self):
        if not isinstance(self.message, str):
            raise ValueError(f"Created a Search Message with something other than string {self.message}")

    def __str__(self):
        return self.message

    @property
    def _sort_order(self):
        return (self.genome_build.name if self.genome_build else "", self.severity, self.message)

    def __lt__(self, other: 'SearchMessage'):
        return self._sort_order < other._sort_order


@dataclass(frozen=True)
class SearchResultGenomeBuildMessages:
    genome_build: GenomeBuild
    messages: List[SearchMessage]

    @cached_property
    def message_summary(self):
        if self.messages:
            return f"In {self.genome_build}<br/>" + "<br/>".join(message.message for message in self.messages)

    @cached_property
    def severity_bs(self):
        if self.messages:
            max_log = max(*(m.severity for m in self.messages), key=lambda s: log_level_to_int(s) )
            return log_level_to_bootstrap(max_log)
        return ""


@dataclass
class SearchResult:
    preview: PreviewData
    messages: List[SearchMessage] = field(default_factory=list)
    sub_name: Optional[str] = None
    match_strength: SearchResultMatchStrength = None
    ignore_genome_build_mismatch: bool = False
    # this is provided automatically by the search framework

    parent: Optional['SearchResponse'] = None
    original_order: Optional[int] = None
    # the above two should be populated by outside code

    @property
    def _sort_order(self):
        return self.genome_build_mismatch, self.original_order

    def __lt__(self, other):
        return self._sort_order < other._sort_order

    @property
    def effective_match_strength(self) -> SearchResultMatchStrength:
        return self.match_strength or SearchResultMatchStrength.STRONG_MATCH

    @cached_property
    def genome_build_relevant_messages(self) -> [SearchMessage]:
        preferred_gb = self.parent.search_input.genome_build_preferred
        if not self.genome_builds:
            return self.messages
        elif preferred_gb not in self.genome_builds:
            return self.messages
        else:
            return [message for message in self.messages if message.genome_build is None or message.genome_build == preferred_gb]

    @property
    def genome_builds_with_messages(self) -> List[SearchResultGenomeBuildMessages]:
        all_results = []
        for genome_build in sorted(self.genome_builds):
            messages = [message for message in self.messages if message.genome_build == genome_build]
            all_results.append(SearchResultGenomeBuildMessages(genome_build=genome_build, messages=messages))
        return all_results

    @cached_property
    def genome_build_mismatch(self) -> bool:
        if not self.ignore_genome_build_mismatch and self.genome_builds and self.parent.search_input.genome_build_preferred not in self.genome_builds:
            return True
        return False

    @property
    def is_perfectly_valid(self):
        return not self.genome_build_relevant_messages and not self.genome_build_mismatch

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
        return self._preview_icon_severity('success')
        # if self.messages:
        #     return self._preview_icon_severity('warning')
        # else:
        #     return self._preview_icon_severity('success')

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


def _default_make_search_result(obj) -> Optional['SearchResult']:
    if isinstance(obj, SearchResult):
        return obj
    else:
        return SearchResult(preview=obj.preview)


@dataclass
class _SearchResultFactory:
    source: Any
    converter: Callable[[Any], SearchResult]
    extra_messages: List[SearchMessage]

    @staticmethod
    def convert(output: Any):
        source: Any
        factory: Callable[[Any], SearchResult] = _default_make_search_result
        extra_messages: List[SearchMessage] = []
        if isinstance(output, (tuple, list)):
            source = output[0]
            for extra in output[1:]:
                if isinstance(extra, str):
                    extra_messages.append(SearchMessage(message=extra))
                elif isinstance(extra, SearchMessage):
                    extra_messages.append(extra)
                elif isinstance(extra, list):
                    for extra_item in extra:
                        if isinstance(extra_item, str):
                            extra_messages.append(SearchMessage(message=extra_item))
                        elif isinstance(extra_item, SearchMessage):
                            extra_messages.append(extra_item)
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


@dataclass()
class SearchResponse:
    search_input: SearchInput
    search_type: PreviewCoordinator
    matched_pattern: bool = False
    sub_name: Optional[str] = None
    admin_only: bool = False
    example: Optional[SearchExample] = None
    results: List[SearchResult] = field(default_factory=list)
    messages: List[SearchMessageOverall] = field(default_factory=list)
    total_count: int = 0

    def __post_init__(self):
        for i, result in enumerate(self.results):
            result.original_order = i
            result.parent = self
        self.results = list(sorted(self.results))

    @property
    def search_messages_with_type(self) -> List[SearchMessagesOverallForType]:
        return [SearchMessagesOverallForType(search_message=message, search_info=self) for message in self.messages]

    @property
    def preview_category(self):
        # django templates don't handle {{ UserPreview.preview_category }} for some unknown reason, works for all models classes
        return self.search_type.preview_category()

    @property
    def preview_icon(self):
        return self.search_type.preview_icon()

    def __lt__(self, other):
        return (self.preview_category, self.sub_name or "") < (other.preview_category, other.sub_name or "")


# consider moving this into variant_search... somehow - maybe provide a transformer in the decorator
def _convert_variant_search_response_to_allele_search_response(variant_response: SearchResponse) -> SearchResponse:
    """
    We perform variant searches, but environments like to consolidate 37 and 38 variants to their alleles.
    This method converts a SearchResponse for Variants into a SearchResponse for Alleles.
    It will get cranky if there are variants without Alleles
    """
    allele_to_variants: Dict[Allele, List[SearchResult]] = defaultdict(list)
    no_allele_variants: List[SearchResult] = []
    allele_results: List[SearchResult] = []
    non_variant_data: List[SearchResult] = []  # probably ManulCreateVariants etc., just call them Allele data

    for result in variant_response.results:
        if obj := result.preview.obj:
            if isinstance(obj, Variant):
                if allele := obj.allele:
                    allele_to_variants[allele].append(result)
                else:
                    no_allele_variants.append(result)
            else:
                non_variant_data.append(result)
        else:
            non_variant_data.append(result)

    for allele, variant_results in allele_to_variants.items():
        genome_builds = set().union(
            *[vr.preview.genome_builds for vr in variant_results if vr.preview.genome_builds])

        all_messages = list(sorted(set().union(*(set(variant.messages) for variant in variant_results if variant.messages))))

        strongest_match = max(variant.match_strength for variant in variant_results)

        ignore_genome_build_mismatch = any(variant.ignore_genome_build_mismatch for variant in variant_results)

        allele_preview = allele.preview
        allele_preview.genome_builds = genome_builds

        allele_results.append(
            SearchResult(
                preview=allele_preview,
                messages=all_messages,
                sub_name=variant_response.sub_name,
                match_strength=strongest_match,
                ignore_genome_build_mismatch=ignore_genome_build_mismatch
            )
        )

    all_results = allele_results + non_variant_data

    top_level_messages = variant_response.messages or []
    if no_allele_variants:
        for no_allele in no_allele_variants:
            no_allele.preview.category = "Allele"
            no_allele.messages.append(SearchMessage("This variant is not linked to an allele", severity=LogLevel.INFO))
            all_results.append(no_allele)

    return SearchResponse(
        search_input=variant_response.search_input,
        search_type=Allele,
        matched_pattern=variant_response.matched_pattern,
        messages=top_level_messages,
        example=variant_response.example,
        sub_name=variant_response.sub_name,
        admin_only=variant_response.admin_only,
        results=all_results,
        total_count=len(all_results)
    )


@dataclass
class SearchCount:
    category: str
    resolved: int = 0
    total: int = 0

    @property
    def excluded_results(self) -> bool:
        return self.total > self.resolved


class SearchResponsesCombined:

    def __init__(self, search_input: SearchInput, responses: List[SearchResponse]):
        self.search_input = search_input
        self.responses = list(sorted(responses))

    @cached_property
    def search_counts(self) -> List[SearchCount]:
        counts: Dict[str, SearchCount] = {}
        for response in self.responses:
            # FIXME if the responses expected return type doesn't match a SearchResult's actual category
            # these counts all go to bunk, probably best to disallow that
            category = response.preview_category
            sc = counts.get(category)
            if sc is None:
                sc = SearchCount(category=category)
                counts[category] = sc
            sc.resolved += len(response.results)
            sc.total += response.total_count

        return list(sorted((sc for sc in counts.values() if sc.resolved), key=lambda x: x.category))

    @property
    def has_excluded_records(self):
        return any(sc for sc in self.search_counts if sc.excluded_results)

    @cached_property
    def results(self) -> List[SearchResult]:
        return list(itertools.chain.from_iterable([response.results for response in self.responses]))

    @cached_property
    def messages(self) -> List[SearchMessagesOverallForType]:
        return list(itertools.chain.from_iterable([response.search_messages_with_type for response in self.responses]))

    def single_preferred_result(self):
        if first(message for message in self.messages if message.severity == LogLevel.ERROR):
            return None

        # If we have ID Matches, ignore other matches for auto-jumping to results
        if id_matches := take(2, (result for result in self.results if result.match_strength == SearchResultMatchStrength.ID_MATCH)):
            if len(id_matches) == 2:
                return None
            if id_matches[0].is_perfectly_valid:
                return id_matches[0]

        # alternatively if we only have 1 match, return it
        elif len(self.results) == 1:
            first_result = self.results[0]
            if first_result.is_perfectly_valid:
                return first_result

        # FIXME support for variant results across genome builds that haven't been merged into alleles

    @property
    def summary(self) -> str:
        return f"{self.search_input.search_string}' calculated {len(self.results)} results."

    @cached_property
    def grouped_search_infos(self) -> List[SearchResponse]:
        return sorted(self.responses, key=lambda sr: (sr.preview_category, sr.sub_name.upper() if sr.sub_name else ""))


def search_data(user: User, search_string: str, classify: bool = False) -> SearchResponsesCombined:
    user_settings = UserSettings.get_for_user(user)
    search_input = SearchInput(
        user=user,
        search_string=clean_string(search_string),
        genome_build_preferred=user_settings.default_genome_build,
        classify=classify
    )
    responses = search_input.search()
    return SearchResponsesCombined(search_input, responses)


def search_receiver(
        search_type: Optional[PreviewCoordinator],
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

            overall_messages: Set[SearchMessageOverall] = set()

            matched_pattern = False
            messages = []
            results = []
            total_count = 0

            try:
                if match := pattern.search(search_input.search_string):
                    matched_pattern = True
                    # as Variants get merged into Alleles, we want to avoid limiting them (except under extreme conditions)
                    limit = 1000 if search_type.preview_category() == "Variant" else 50
                    for result in func(SearchInputInstance(expected_type=search_type, search_input=search_input, match=match)):
                        if result is None:
                            raise ValueError(f"Search {sender.__name__} returned None")
                        elif isinstance(result, SearchMessageOverall):
                            overall_messages.add(result)
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

                                    results.append(search_result)
                                    limit -= 1

            except Exception as e:
                # TODO, determine if the Exception type is valid for users or not
                overall_messages.append(SearchMessageOverall(str(e), severity=LogLevel.ERROR))
                print(f"Error handling search_receiver on {func}")
                report_exc_info()

            response = SearchResponse(
                search_input=search_input,
                search_type=search_type,
                admin_only=admin_only,
                sub_name=sub_name,
                example=example,
                matched_pattern=matched_pattern,
                results=results,
                messages=list(sorted(overall_messages)),
                total_count=total_count
            )

            if settings.PREFER_ALLELE_LINKS and response.search_type.preview_category() == "Variant":
                response = _convert_variant_search_response_to_allele_search_response(response)

            return response

        search_signal.connect(search_func)
        return search_func

    return _decorator
