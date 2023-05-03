import operator
import re
import itertools
from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from functools import cached_property, reduce
from re import IGNORECASE
from typing import List, Set, Optional, Type, Pattern, Callable, Any, Match, Union, Dict, Iterable

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet, Q
from django.dispatch import Signal
from more_itertools import take

from library.enums.log_level import LogLevel
from library.log_utils import report_exc_info, report_message, log_level_to_int, log_level_to_bootstrap
from library.preview_request import PreviewCoordinator, PreviewData, PreviewModelMixin
from library.utils import clean_string, first
from snpdb.models import UserSettings, GenomeBuild, Variant, Allele


search_signal = Signal()
HAS_ALPHA_PATTERN = re.compile(r"[a-zA-Z]")
HAS_ANYTHING = re.compile(r".")
HAS_3_ALPHA_MIN = re.compile(r"[a-zA-Z]{3,}")
HAS_3_ANY = re.compile(r"\S{3,}")
_SPLIT_GAPS = re.compile(r"[\s,]+")

MAX_VARIANT_RESULTS = 100
MAX_RESULTS_PER_SEARCH = 50


"""
Important: If you just want to create a new search, jump to the @search_receiver and read the documention there

Key components of search

SearchInput: Handles all the adjustable parameters of search, the user who iniated the search, their preferred genome build, the text they're searching for etc.
SearchInputInstance: As above but stores the regex match for a specific search (useful if a search's required pattern has capture groups)
@search_receiver: A decoration on a method that takes a SearchInputInstance and gives a SearchResponse. There's a bit of decorator magic to do the conversion.
SearchMessageOverall: A message that's against a specific search, but not against a specific search result
SearchMessage: A message that's against a specific search result
SearchResult: Relates to one row in the search result listing, keeps information such as preview data (icon, category, labels), and also gets some meta information from the SearchResponse that created it
SearchResponse: A combination of SearchResults and SearchMessageOveralls
SearchResponseCombined: The complete output of a search across all teh different types
"""


@dataclass(frozen=True)
class SearchInput:
    user: User
    search_string: str
    genome_build_preferred: GenomeBuild
    classify: bool = False
    """
    Is this coming from the classify by c.HGVS form
    """

    def matches_pattern(self, pattern: Union[str, Pattern]) -> Match:
        if isinstance(pattern, str):
            return re.match(pattern, self.search_string, IGNORECASE)
        else:
            return pattern.match(self.search_string)

    @property
    def search_words(self) -> List[str]:
        """
        Split on white space and commas in the search input - use when filtering on name
        :return:
        """
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

        # TODO, this would be a good place to handle roman numerals and arabic as the inverse (per Ontology search)
        # have a parameter to say if the substitution is desired
        if words := self.search_words:
            qs: List[Q] = []
            for word in words:
                qs.append(Q(**{f"{field_name}__{test}": word}))
            return reduce(operator.and_, qs)
        else:
            raise ValueError("No tokens found in search, can't generate q_words")

    def search(self) -> List['SearchResponse']:
        """
        Execute the search by firing off the search_signal (which is connected to via @search_receiver)
        :return:
        """
        valid_responses: List[SearchResponse] = []
        response_tuples = search_signal.send_robust(sender=SearchInput, search_input=self)
        for caller, response in response_tuples:
            if response:
                if response is None:
                    # we may occasionally skip searches for some reasons, such as we're here from classify
                    # so we only care about variant matches
                    continue
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
        """
        Provide the list of genome builds that should be included in the search
        """
        return set(GenomeBuild.builds_with_annotation().all())


@dataclass(frozen=True)
class SearchInputInstance:
    """
    When invoking individual @search_receivers, customise the SearchInputInstance to provide the pattern match.
    Useful for @search_receivers that have a capture group, not real bonus otherwise
    """

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

    @property
    def genome_build_preferred(self) -> GenomeBuild:
        return self.search_input.genome_build_preferred

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
    search_info: Optional['SearchResponse'] = None

    @property
    def _sort_order(self):
        return log_level_to_int(self.severity), self.message

    def __lt__(self, other) -> bool:
        return self._sort_order < other._sort_order

    def __post_init__(self):
        if not isinstance(self.message, str):
            raise ValueError(f"Created a Search Message with something other than string {self.message}")

    def with_response(self, search_response: 'SearchResponse'):
        return SearchMessageOverall(message=self.message, severity=self.severity, search_info=search_response)


class SearchResultMatchStrength(int, Enum):
    """
    Each @search_receiver will have a default search strength.
    Just matching on parts of the name should return a fuzzy_match
    Matching using logic (such as on a variant HGVS) is a strong match
    Straight up providing the ID, e.g. PMID:3453453 or BRCA2 is an ID match
    Helps determine if we should auto jump to any result
    """
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
        return self.genome_build.name if self.genome_build else "", log_level_to_int(self.severity), self.message

    def __lt__(self, other: 'SearchMessage'):
        return self._sort_order < other._sort_order

    def with_genome_build(self, genome_build: GenomeBuild) -> 'SearchMessage':
        return SearchMessage(message=self.message, severity=self.severity, genome_build=genome_build)

    @staticmethod
    def highest_severity(iterable: Iterable['SearchMessage']) -> LogLevel:
        return max((m.severity for m in iterable), key=lambda sm: log_level_to_int(sm), default=LogLevel.INFO)


@dataclass(frozen=True)
class SearchResultGenomeBuildMessages:
    """
    Summary of genome builds in the response, and messages for that genome build.
    Note that genome build agnostic messages will return an empty array of these
    """
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
    """
    Represents a single record found in a search
    """

    preview: PreviewData
    """
    The actual object returned by the search
    """

    messages: List[SearchMessage] = field(default_factory=list)
    """
    Any validation messages associated with the preview
    """

    match_strength: SearchResultMatchStrength = None
    """
    How strong was the match, will default to the @search_receiver's max strength
    """

    ignore_genome_build_mismatch: bool = False
    """
    Very specific usage, if searching on g. where the genome build is built into the search data,
    don't warn the user that the value returned isn't their default genome build
    """

    parent: Optional['SearchResponse'] = None
    """
    Populated automatically, don't set this in @search_receivers
    """

    original_order: Optional[int] = None
    """
    Populated automatically, don't set this in @search_receivers.
    Allows the SearchResponse to maintain the order of results as provided to it, but provide higher level sorting
    based on preferred genome build or other factors
    """

    @property
    def sub_name(self) -> str:
        if self.parent:
            return self.parent.sub_name

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
        """
        What validation messages are relevant to show for this record, e.g. if it's a variant found using both 37 & 38
        and the user's preferred genome build is 38, this will return messages associated to no genome build and to
        38 (ignoring any associated to 37)
        """
        preferred_gb = self.parent.search_input.genome_build_preferred
        if not self.genome_builds:
            return self.messages
        elif preferred_gb not in self.genome_builds:
            return self.messages
        else:
            return [message for message in self.messages if message.genome_build is None or message.genome_build == preferred_gb]

    @property
    def genome_builds_with_messages(self) -> List[SearchResultGenomeBuildMessages]:
        """
        Lists genome build with associated messages
        """
        all_results = []
        for genome_build in sorted(self.genome_builds):
            messages = [message for message in self.messages if message.genome_build == genome_build]
            all_results.append(SearchResultGenomeBuildMessages(genome_build=genome_build, messages=messages))
        return all_results

    @cached_property
    def genome_build_mismatch(self) -> bool:
        """
        Returns if the result is not genome build agnostic, is not associated with the user's preferred genome build, and the result doesn't
        ignore genome build mismatches
        """
        if not self.ignore_genome_build_mismatch and self.genome_builds and self.parent.search_input.genome_build_preferred not in self.genome_builds:
            return True
        return False

    @property
    def is_perfectly_valid(self):
        """
        Make sure there are no relevant errors or genome build mismatch or if this is an operation that should probably not occur instantly.
        Used to determine if a result can be auto-redirected to
        """
        return SearchMessage.highest_severity(self.genome_build_relevant_messages) not in ('W', 'E') and not self.genome_build_mismatch and not self.preview.is_operation

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
        # coloring the icon became too busy, we're already showing info/warning/error icons against the record
        # no need to re-color the icon
        return self._preview_icon_severity('success')
        # if self.messages:
        #     return self._preview_icon_severity('warning')
        # else:
        #     return self._preview_icon_severity('success')

    @property
    def search_type(self) -> str:
        return self.preview.category

    @property
    def genome_builds(self) -> Optional[Set[GenomeBuild]]:
        return self.preview.genome_builds

    @property
    def annotation_consortia(self) -> Optional[List['AnnotationConsortia']]:
        if consortia := self.preview.annotation_consortia:
            return list(sorted(consortia))


def _default_make_search_result(obj) -> Optional['SearchResult']:
    """
    This is the default converter for making an object into a SearchResult (based on what is yielded by a
    @search_receiver. Handling of iterables, and tuple data is handled outside of this.
    Note that a @search_receiver can return a factory Callable to replace this.
    """
    if isinstance(obj, SearchResult):
        return obj
    elif isinstance(obj, PreviewData):
        return SearchResult(preview=obj)
    else:
        return SearchResult(preview=obj.preview)


@dataclass
class _SearchResultFactory:
    """
    Is in charge of converting raw output from a @search_receiver into SearchResults
    """

    source: Any
    """
    Could be a
    * QuerySet, a PreviewData, an instance of a PreviewModelMixin or list
    * A tuple where the first item is any of the above, and the subsequent items are:
    ** str, List[str], SearchMessage, List[SearchMessage] - to be converted to SearchMessage and applies to every
    value from the first index.
    ** A callable that takes a value from the first item of the tuple and creates a SearchResult (handy for providing
    extra validation).
    """

    converter: Callable[[Any], SearchResult]
    """
    A callable that takes a value from the first item of the tuple and creates a SearchResult (handy for providing
    extra validation)
    """

    extra_messages: List[SearchMessage]
    """
    If str, List[str] etc is provided in a tuple yielded from a @search_receiver, the messages apply to every object
    that came from the first parameter
    """

    @staticmethod
    def convert(output: Any):
        source: Any
        factory: Callable[[Any], SearchResult] = _default_make_search_result
        extra_messages: List[SearchMessage] = []

        if isinstance(output, tuple):
            # if we got a tuple or a list, the first entry is the "source" and the subsequent values
            # will be messages or a converter
            source = output[0]
            for extra in output[1:]:
                # loop through the
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
                elif isinstance(extra, Callable):
                    factory = extra
                else:
                    raise ValueError(f"Search method yielded extra {extra}, expected str, list[str] or callable")
        else:
            # should be a QuerySet (of PreviewModelMixin), list (of PreviewModelMixin instances), or a single instance of one
            source = output

        return _SearchResultFactory(source=source, converter=factory, extra_messages=extra_messages)

    @property
    def _iterable(self):
        """
        Turn source into a splicable iteratable object
        :return:
        """
        if isinstance(self.source, (list, QuerySet)):
            return self.source
        else:
            return [self.source]

    def __len__(self) -> int:
        """
        Returns the maximum number of results this factory could return
        (importantly returns .count() on QuerySets). Allows us to return a max of 50 results, while also reporting
        the maximum number of results that were returned.
        """
        if isinstance(self.source, list):
            return len(self.source)
        elif isinstance(self.source, QuerySet):
            return self.source.count()
        else:
            return 1

    def iterate(self, limit: int = MAX_RESULTS_PER_SEARCH) -> Iterable[SearchResult]:
        """
        start spitting out SearchResponse up to a limit of records
        :param limit: The maximum number of results to return, if multiple things were yielded keep reducing the limit based on how many
        records were returned by previous _SearchResultFactories
        :return: An iterator that will return 0 to limit SearchResult
        """
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
    """
    For a @search_receiver, provide an example of how the data can be searched. Is shown to the user.
    """

    note: Optional[str] = None
    """
    English explanation of what text will activate this search
    """

    examples: Optional[List[str]] = None
    """
    An example of input that will activate the search
    """


@dataclass()
class SearchResponse:
    """
    The output of a @search_receiver, e.g. all the MONDO terms found by a MONDO name search. Variant search is
    broken up into multiple searches e.g. c.HGVS, Variant Coordinate, and each will return a SearchResponse
    """

    search_input: SearchInput
    """
    What was the input that resulted in this output
    """

    search_type: PreviewCoordinator
    """
    What's our default icon and category
    """

    matched_pattern: bool = False
    """
    Did the search input match the basic pattern required for this search to fire. If false, there should be no
    results or messages_overall. If there was no search input we still get all the search_responses
    """

    sub_name: Optional[str] = None
    """
    The sub_name as taken from the @search_receiver
    """

    admin_only: bool = False
    """
    Is this search for admin users only
    """

    example: Optional[SearchExample] = None
    """
    The example to show to the user
    """

    results: List[SearchResult] = field(default_factory=list)
    """
    Actual search results if matched_pattern
    """

    messages_overall: List[SearchMessageOverall] = field(default_factory=list)
    """
    Messages about the search overall that aren't linked to a single record.
    Handy if the input is broken
    """

    total_count: int = 0
    """
    How many search results were calculated (should be equal to or larger than the length of results)
    """

    def __post_init__(self):
        # once you create a SearchResponse, don't mutate it as __post_init__ fixes a few things
        # can't make it frozen due to
        for i, result in enumerate(self.results):
            result.original_order = i
            result.parent = self
            result.messages = list(sorted(result.messages))
        self.results = list(sorted(self.results))
        self.messages_overall = [om.with_response(self) for om in self.messages_overall]

    @property
    def preview_category(self) -> str:
        # overall category (the category on individual results within a SearchResponse can be of different categories
        # but this is not recommended
        return self.search_type.preview_category()

    @property
    def preview_icon(self) -> str:
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
    non_variant_data: List[SearchResult] = []  # probably ManualCreateVariants etc

    for result in variant_response.results:
        if obj := result.preview.obj:
            if isinstance(obj, Variant):
                if allele := obj.allele:
                    allele_to_variants[allele].append(result)
                else:
                    # hack variant icon to be Allele icon even if there's no allele
                    result.preview_icon = Allele.preview_icon()
                    no_allele_variants.append(result)
            else:
                non_variant_data.append(result)
        else:
            non_variant_data.append(result)

    for allele, variant_results in allele_to_variants.items():
        genome_builds = set().union(
            *[vr.preview.genome_builds for vr in variant_results if vr.preview.genome_builds])

        all_messages = []
        for variant in variant_results:
            all_messages += [r.with_genome_build(first(variant.genome_builds)) for r in variant.messages]

        all_messages = list(sorted(set(all_messages)))

        strongest_match = max(variant.match_strength for variant in variant_results)

        ignore_genome_build_mismatch = any(variant.ignore_genome_build_mismatch for variant in variant_results)

        allele_preview = allele.preview
        allele_preview.genome_builds = genome_builds

        allele_results.append(
            SearchResult(
                preview=allele_preview,
                messages=all_messages,
                match_strength=strongest_match,
                ignore_genome_build_mismatch=ignore_genome_build_mismatch
            )
        )

    all_results = allele_results + non_variant_data

    top_level_messages = variant_response.messages_overall or []
    if no_allele_variants:
        for no_allele in no_allele_variants:
            no_allele.preview.category = "Allele"
            no_allele.messages.append(SearchMessage("This variant is not yet linked to an allele", severity=LogLevel.INFO))
            all_results.append(no_allele)

    return SearchResponse(
        search_input=variant_response.search_input,
        search_type=Allele,
        matched_pattern=variant_response.matched_pattern,
        messages_overall=top_level_messages,
        example=variant_response.example,
        sub_name=variant_response.sub_name,
        admin_only=variant_response.admin_only,
        results=all_results,
        total_count=len(all_results)
    )


@dataclass
class SearchCount:
    """
    Used to provide a search result count summary in the nav bar
    """

    category: str
    resolved: int = 0
    total: int = 0

    @property
    def excluded_results(self) -> bool:
        return self.total > self.resolved


class SearchResponsesCombined:
    """
    Is the complete output of a SearchInput across all searches
    """

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
    def messages_overall(self) -> List[SearchMessageOverall]:
        return list(itertools.chain.from_iterable([response.messages_overall for response in self.responses]))

    def single_preferred_result(self):
        if first(message for message in self.messages_overall if log_level_to_int(message.severity) >= log_level_to_int(LogLevel.WARNING)):
            return None

        # If we have ID Matches, ignore other matches for auto-jumping to results
        if id_matches := take(2, (result for result in self.results if result.match_strength == SearchResultMatchStrength.ID_MATCH)):
            if len(id_matches) == 2:
                return None
            if id_matches[0].is_perfectly_valid:
                return id_matches[0]
            return None

        correct_genome_build_results = [result for result in self.results if not result.genome_build_mismatch]

        # alternatively if we only have 1 match, return it
        if len(correct_genome_build_results) == 1:
            first_result = correct_genome_build_results[0]
            if first_result.is_perfectly_valid and first_result.match_strength >= SearchResultMatchStrength.STRONG_MATCH:
                return first_result

    @property
    def summary(self) -> str:
        # is recorded in the EventLog
        return f"{self.search_input.search_string}' calculated {len(self.results)} results."


def search_data(user: User, search_string: str, classify: bool = False) -> SearchResponsesCombined:
    """
    Call externally to do all the magic of searching
    :param user: The user performing the search
    :param search_string: The text that we entered
    :param classify: Have we come from a user wanting to classify a given HGVS
    :return: SearchResponseCombined
    """
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
    Wrap around a Callable[[SearchInputInstance], Generator[Any]]
    The wrapped method can return:
    A single PreviewData, Model[PreviewDataMixin], SearchResult
    A QuerySet or List of PreviewDataMixin
    A tuple where the first element is one of the above, and the subsequent items:
    * SearchMessage or str that will be turned into SearchMessage and/or
    * A Callable[[T], SearchResult] where T is what the first element returns (or the first element is a QuerySet of).

    For some examples:

    def super_simple_search(search_input: SearchInputInstance):
        yield Model.objects.filter(xxx)

    def simple_search(search_input: SearchInputInstance):
        for x in filtered_results:
            yield x

    def medium_search(search_input: SearchInputInstance):
        for x in filtered_results:
            messages: List[str] = validate(x)
            yield x, messages

    def complicated_search(search_input: SearchInputInstance):
        def validate_data(data):
            ...
            return SearchResult(preview_of_data(data), validation_message_for_data)
        return Data.objects.filter(xxx), validate_data

    The complicated_search example doesn't iterate over every result, but provide a method that will annotate a result with
    validation messages. This is useful for lazily validating only relevant results when we're going to limit the results
    returned to 50 for example.

    Note that after the decoration, the wrapped method signature will take a SearchInput and return a SearchResponse object.

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
            results = []
            total_count = 0

            if search_input.classify and search_type.preview_category() != "Variant":
                return None

            try:
                if match := pattern.search(search_input.search_string):
                    matched_pattern = True
                    # as Variants get merged into Alleles, we want to avoid limiting them (except under extreme conditions)
                    limit = MAX_VARIANT_RESULTS if search_type.preview_category() == "Variant" else MAX_RESULTS_PER_SEARCH
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

                                    results.append(search_result)
                                    limit -= 1

            except Exception as e:
                # TODO, determine if the Exception type is valid for users or not
                overall_messages.add(SearchMessageOverall(str(e), severity=LogLevel.ERROR))
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
                messages_overall=list(sorted(overall_messages)),
                total_count=total_count
            )

            if settings.PREFER_ALLELE_LINKS and response.search_type.preview_category() == "Variant":
                try:
                    response = _convert_variant_search_response_to_allele_search_response(response)
                except:
                    report_exc_info()
                    response.messages_overall.append(SearchMessageOverall("Unexpected error when attempting to convert Variant results into Allele results"))

            return response

        search_signal.connect(search_func)
        return search_func

    return _decorator
