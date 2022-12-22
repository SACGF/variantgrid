import datetime
import itertools
import re
import typing
from dataclasses import dataclass, field
from datetime import timedelta
from enum import Enum
from typing import Optional, Iterable, List, Dict, Set, Union, Any, Iterator, Tuple

from Bio import Entrez, Medline
from django.db import models
from django.db.models import TextField
from django.urls import reverse
from django.utils.timezone import now
from django_extensions.db.models import TimeStampedModel

from library.log_utils import report_exc_info
from library.utils import JsonObjType, first

"""
Has the model for Citation as well as the methods for populating them.
Short version

CitationFetchRequest.fetch_all_now(["PMID:1231232", "PMID:2342323"]).citations
# returns an array of populated citations

Going via CitationFetchRequest should generally be the only way you interact with Citations
"""


DATE_PUBLISHED_RE = re.compile("^([0-9]+).*?$")


def get_year_from_date(date_published: str) -> str:
    if not date_published:
        return ""
    year = ""
    if year_match := DATE_PUBLISHED_RE.fullmatch(date_published):
        year = year_match.group(1)
    return year


class CitationSource(models.TextChoices):
    PUBMED = 'PMID', 'PubMed'
    NCBI_BOOKSHELF = 'Bookshelf ID', 'NCBI Bookshelf'
    PUBMED_CENTRAL = 'PMCID', 'PubMedCentral'

    @staticmethod
    def from_legacy_code(code: str) -> Optional['CitationSource']:
        """
        Turns any kind of possible citation prefix into a CitationSource
        :param code: String representing PubMed, PubMedCentral or NCBI's Bookshelf
        :return: A CitationSource if valid, None otherwise
        """
        return {
            "P": CitationSource.PUBMED,
            "PMID": CitationSource.PUBMED,
            "PUBMED": CitationSource.PUBMED,

            "C": CitationSource.PUBMED_CENTRAL,
            "PMCID": CitationSource.PUBMED_CENTRAL,
            "PubMedCentral": CitationSource.PUBMED_CENTRAL,

            "N": CitationSource.NCBI_BOOKSHELF,
            "NBK": CitationSource.NCBI_BOOKSHELF,
            "NCBIBOOKSHELF": CitationSource.NCBI_BOOKSHELF,
            "BOOKSHELF": CitationSource.NCBI_BOOKSHELF,
            "BOOKSHELF ID": CitationSource.NCBI_BOOKSHELF
        }.get(code.upper())


class EntrezDbType(str, Enum):
    PUBMED = "pubmed"
    PUBMED_CENTRAL = "pmc"
    BOOKSHELF = "books"


class Citation(TimeStampedModel):
    id = TextField(primary_key=True)
    """
    3rd party's unique ID for the citation, e.g. PMID:234434, BookShelf ID:NBK52333
    """

    old_id = models.IntegerField(null=True, blank=True)
    """
    Integer ID when Citaiton used that as its primary key.
    Only needed for Classifications that haven't been re-validated, or historic ClassificationModifications
    """

    source = models.CharField(max_length=16, choices=CitationSource.choices)
    """
    Which kind of citation, e.g. PubMed, PubMedCentral
    """

    index = models.TextField(null=True)
    """
    The non prefix : part of the citation, often contains redundant of its own e.g.
    PMCID:PMC42343 index will be "PMC42343"
    """

    data_json = models.JSONField(null=True, blank=True)
    """
    The data retrieved from Entrez or similar
    """

    last_loaded = models.DateTimeField(null=True, blank=True)
    """
    When this citation data was last requested from Entrez (if it has been).
    Note if this migration was copied over from the old model, this will be null.
    If last_loaded is null, title, journal, journal_short etc should all be blank.
    """

    error = models.TextField(null=True, blank=True)
    """
    Error (if any) when last attempted to retrieve this record
    """

    # details about the citation directly
    title = models.TextField(null=True, blank=True)
    journal = models.TextField(null=True, blank=True)
    journal_short = models.TextField(null=True, blank=True)
    year = models.TextField(null=True, blank=True)
    authors = models.TextField(null=True, blank=True)
    authors_short = models.TextField(null=True, blank=True)
    abstract = models.TextField(null=True, blank=True)

    @property
    def id_pretty(self):
        return self.id.replace(":", ": ")

    @property
    def _first_and_single_author(self) -> Tuple[str, str]:
        first_author = self.authors_short
        single_author = self.authors_short == self.authors
        if self.authors and not first_author:
            author_list = self.authors.split(',')
            if author_list:
                first_author = author_list[0]
            single_author = len(author_list) == 1
        return first_author, single_author

    @property
    def first_author(self):
        return self._first_and_single_author[0]

    @property
    def single_author(self):
        return self._first_and_single_author[1]

    @property
    def title_safe(self):
        if self.error:
            return "Could not retrieve citation"
        elif not self.title:
            return "Have not retrieved citation"
        else:
            return self.title

    def __str__(self):
        return self.id

    @property
    def get_external_url(self):
        if self.source == CitationSource.PUBMED:
            return f"https://www.ncbi.nlm.nih.gov/pubmed/{self.index}"
        elif self.source == CitationSource.PUBMED_CENTRAL:
            return f"https://www.ncbi.nlm.nih.gov/pmc/?term={self.index}"
        else:
            return f"https://www.ncbi.nlm.nih.gov/books/{self.index}"

    def get_absolute_url(self):
        return reverse('view_citation', kwargs={"citation_id": self.id})

    def blank_out(self):
        """
        Removes all the derived data in a Citation
        Generally call before populating data
        """
        self.title = None
        self.journal = None
        self.journal_short = None
        self.year = None
        self.authors = None
        self.authors_short = None
        self.abstract = None


@dataclass(frozen=True)
class CitationIdNormalized:
    """
    A class that represents a validated Citation.
    A CitationID should always be valid in format, though it may reference IDs that don't exist.
    Use a static factory method such as CitationIdNormalized.from_parts and CitationIdNormalized.normalize_id etc
    so validation is performed.
    """

    source: CitationSource
    index: str

    @property
    def full_id(self) -> str:
        """
        Produce the ID as it would be saved in our database and referred to by 3rd parties
        :return: A string identifying the citation
        """
        if self.source == CitationSource.PUBMED:
            return f"PMID:{self.index}"
        elif self.source == CitationSource.PUBMED_CENTRAL:
            return f"PMCID:{self.index}"
        elif self.source == CitationSource.NCBI_BOOKSHELF:
            return f"Bookshelf ID:{self.index}"
        else:
            raise ValueError(f"Unexpected citation source {self.source}")

    CITATION_SPLIT_RE = re.compile(r"(?P<prefix>[a-z ]+)\s*(?P<semicolon>:)?\s*(?P<number_prefix>[a-z]+)?(?P<number>[0-9]+)", re.IGNORECASE)
    """
    Used to split citation strings into parts, e.g. PMCID: PMC32432232 prefix=PMCID, semicolon=: number_prefix=PMC, number=32432232
    """

    NUMER_STRIP_RE = re.compile(r"(?P<number_prefix>[a-z]+)?(?P<number>[0-9]+)", re.IGNORECASE)
    """
    Extract the pure number from index, e.g. PMC32432232 gives you PMC and 32432232
    """

    def __lt__(self, other):
        if self.source < other.source:
            return True
        if self.source == other.source:
            return self.index.rjust(10, '0') < other.index.rjust(10, '0')
        return False

    @staticmethod
    def from_parts(source: Union[str, CitationSource], index: Union[str, int]):
        index = str(index)
        use_source = CitationSource.from_legacy_code(source)
        if not use_source:
            print(f"Unexpected source {source}")
            raise ValueError(f"Unexpected source for Citation ID {source}")

        if match := CitationIdNormalized.NUMER_STRIP_RE.match(index):
            index = match.group('number')
        else:
            raise ValueError(f"Cannot convert {index} to a number")

        if use_source == CitationSource.PUBMED_CENTRAL:
            index = f"PMC{index}"
        elif use_source == CitationSource.NCBI_BOOKSHELF:
            index = f"NBK{index}"
        return CitationIdNormalized(use_source, index)

    @staticmethod
    def from_citation(citation: Citation):
        return CitationIdNormalized.from_parts(source=citation.source, index=citation.index)

    @staticmethod
    def normalize_id(citation_id: Union[str, 'CitationIdNormalized']) -> 'CitationIdNormalized':
        if isinstance(citation_id, CitationIdNormalized):
            return citation_id
        citation_id = citation_id.strip().upper()
        if parts := CitationIdNormalized.CITATION_SPLIT_RE.match(citation_id):
            prefix = parts.group('prefix')
            number = parts.group('number')

            return CitationIdNormalized.from_parts(prefix, number)
        raise ValueError(f"Unable to parse Citation ID {citation_id}")

    def get_or_create(self) -> Citation:
        """
        Gets or creates a Citation for this CitationID
        :return: A Citation
        """
        citation, _ = Citation.objects.get_or_create(
            id=self.full_id,
            defaults={
                "source": self.source,
                "index": self.index
            }
        )
        return citation

    def for_bulk_create(self) -> Citation:
        """
        If using Citation.objects.bulk_create(...., ignore_conflicts=True) use this method
        for objects, populates created, modified which bulk_create wont do
        :return: A Citation for use in bulk_create
        """
        return Citation(
            id=self.full_id,
            source=self.source,
            index=self.index,
            created=now(),
            modified=now()
        )


@dataclass
class CitationFetchEntry:
    """
    Represents the request for a populated Citation.
    Citations end up in 3 states:
    A valid ID, where we don't have a Citation record in the database
    A valid ID, where we do have a Citation record in the database, but it hasn't been populated from Entrez
    A valid ID, where we do have a Citation record in the database and it has been loaded
    An invalid ID
    """

    normalised_id: Optional[CitationIdNormalized] = None
    """
    The requested ID (derived or provided directly)
    """

    citation: Optional[Citation] = None
    """
    The corresponding Citation record
    """

    requested_ids: Set[Any] = field(default_factory=set)
    """
    Used to link what was requested to what we respond with.
    e.g. if an API requested an array of ID strings including "PubMed:   2343223", and the resulting normalised_id, citation don't match
    we can use requested_ids to work out that the resulting "PMID:2343223" matches the requested "PubMed:   2343223"
    """

    fetched: bool = False
    """
    Used internally within CitationFetchRequest to keep track of what records have been processed
    """

    error: Optional[str] = None
    """
    Error if any with the fetch entry
    """

    def record_request_id(self, obj: Any):
        """
        Attempts to record the object used to make the request, must be Hashable to be recorded
        """
        if isinstance(obj, typing.Hashable):
            self.requested_ids.add(obj)

    def should_refresh(self, refresh_before: datetime) -> bool:
        """
        Does this record need to be reloaded from Entrez
        :param refresh_before: If the citation was loaded before this date - or never loaded, reload it
        :return: A boolean indicating if a refresh is required
        """
        if last_loaded := self.citation.last_loaded:
            return last_loaded < refresh_before
        else:
            return True

    def to_json(self) -> JsonObjType:
        if self.citation:
            data = {key: value for key, value in vars(self.citation).items() if isinstance(value, str)}
            data["external_url"] = self.citation.get_external_url
            data["first_author"] = self.citation.first_author
            data["single_author"] = self.citation.single_author
        else:
            data = {"error": self.error}
        data['requested_using'] = [requested_id for requested_id in self.requested_ids if isinstance(requested_id, str)]
        return data


CitationRequest = Union[str, CitationIdNormalized, Citation, int, 'VCDbRefDict', 'DbRefRegexResult']
"""
All the ways a citation can be requested
str : A string (that we will attempt to normalise) e.g. "PMID: 23423423"
CitationIdNormalized: ID is already in a known format
Citation: Citation model, might still need to be populted from Entrez
int: Refers to the old_id in Citation
VCDbRefDict: Expects a dict with a key of 'id' (which then acts the same way as str)
DbRefRegexResult: Expects an attribute 'id_fixed' (which then acts the same way as str)
"""


class CitationFetchResponse:
    """
    Wraps the citations returned from a CitationFetchRequest
    """

    def __init__(self, entries: List[CitationFetchEntry]):
        self.all_entries = entries
        requested_to_fetch: Dict[Any, CitationFetchEntry] = dict()
        for fetch in entries:
            for requested_id in fetch.requested_ids:
                requested_to_fetch[requested_id] = fetch
        self.requested_to_fetch = requested_to_fetch

    @property
    def all_citations(self) -> List[Citation]:
        return [entry.citation for entry in self.all_entries]

    def for_requested(self, citation_id: CitationRequest) -> Citation:
        """
        Return (a hopefully populated) Citation for the ID that was used to request it
        """
        return self.requested_to_fetch.get(citation_id).citation

    def __iter__(self) -> Iterator[CitationRequest]:
        return iter(self.all_entries)

    def __len__(self) -> int:
        return len(self.all_entries)

    @property
    def first_citation(self) -> Citation:
        return first(self.all_entries).citation

    def to_json(self):
        return [entry.to_json() for entry in self.all_entries]


class CitationFetchRequest:
    """
    Use to request a batch of populated Citations. Attempts to do so with the minimum number of database and network calls.
    """

    TOP_LEVEL_NBK_RE = re.compile("^NBK[0-9]+$")
    """
    Used to identify which response from requesting NBK represents the book (and not just a chapter)
    """

    def __init__(self, cache_age: Optional[timedelta] = None):
        self.id_to_fetch: Dict[CitationIdNormalized, CitationFetchEntry] = dict()
        """
        Keeps a dict of all normalized IDs to FetchEntries but only for records
        """

        self.error_fetches: List[CitationFetchEntry] = list()
        """
        FetchEntries that are too invalid to request from Entrez
        """

        self.refresh_before = now() - (cache_age if cache_age is not None else timedelta(days=730))  # reload every 2 years

    @staticmethod
    def fetch_all_now(citation_ids: Iterable[CitationRequest], cache_age: Optional[timedelta] = None) -> CitationFetchResponse:
        """
        Return Citations that we attempt to populate
        :param citation_ids: A series of requests, be they strings, or formatted dicts (note that dicts for non citations will be ignored)
        :param cache_age: Reload any citations older than the given date, is 2 years by default (just in case there's a correction)
        """
        cfr = CitationFetchRequest(cache_age=cache_age)
        for citation_id in citation_ids:
            cfr._queue(citation_id)
        cfr._fetch_queue()
        return CitationFetchResponse(cfr._all_citation_fetches + cfr.error_fetches)

    @staticmethod
    def get_unfetched_citations(citation_ids: Iterable[CitationRequest]) -> List[Citation]:
        """
        Return Citation objects, but don't attempt to populate from Entrez (records might already be populated)
        """
        cfr = CitationFetchRequest()
        for citation_id in citation_ids:
            cfr._queue(citation_id)
        cfr._load_citation_stubs()
        return [cfe.citation for cfe in cfr._all_citation_fetches if cfe.citation]

    def _queue(self, citation_id: CitationRequest) -> CitationFetchEntry:
        citation: Optional[Citation] = None
        normalized: Optional[CitationIdNormalized] = None
        # handle a WHOLE slew of inputs

        try:
            # if just a number, assume it's the old Citation ID which should match a record's old_id value
            # migrate away from this when possible
            if isinstance(citation_id, int) or (isinstance(citation_id, str) and citation_id.isnumeric()):
                citation = Citation.objects.filter(old_id=int(citation_id)).first()
                if not citation:
                    # something's gone wrong, an old_Id was provided, but we have no copy of it
                    invalid_id_fetch = CitationFetchEntry()
                    invalid_id_fetch.record_request_id(citation_id)
                    invalid_id_fetch.error = "Internal Error: could not load citation"
                    self.error_fetches.append(invalid_id_fetch)
                    return invalid_id_fetch

                normalized = CitationIdNormalized.normalize_id(citation.pk)

            # Already a Citation, easy
            elif isinstance(citation_id, Citation):
                citation = citation_id
                normalized = CitationIdNormalized.normalize_id(citation.id)

            # Already a CitationIdNormalized, easy
            elif isinstance(citation_id, CitationIdNormalized):
                normalized = citation_id

            # A string, assume it's a valid citation reference, just might need normalisation
            elif isinstance(citation_id, str):
                normalized = CitationIdNormalized.normalize_id(citation_id)

            # A dict, assume it's a VCDbRefDict where ["id"] will be a valid citation reference
            elif isinstance(citation_id, dict):
                normalized = CitationIdNormalized.normalize_id(citation_id.get("id"))

            # A DbRefRegexResult,
            elif hasattr(citation_id, 'id_fixed'):
                normalized = CitationIdNormalized.normalize_id(citation_id.id_fixed)

        except ValueError:
            # Most likely asked for db_ref not for citations... don't add these to results
            # Invalid ID, don't even attempt to fetch this
            # asking the FetchEntry if is_valid will return False
            invalid_id_fetch = CitationFetchEntry()
            invalid_id_fetch.record_request_id(citation_id)
            invalid_id_fetch.error = f"{citation_id} could not load, invalid ID"
            return invalid_id_fetch

        existing = self.id_to_fetch.get(normalized)
        entry: CitationFetchEntry
        if existing:
            entry = existing
        else:
            entry = CitationFetchEntry(
                normalised_id=normalized,
                citation=citation
            )
        entry.record_request_id(citation_id)

        self.id_to_fetch[normalized] = entry
        return entry

    def _fetch_now(self, citation_id: CitationRequest) -> CitationFetchEntry:
        entry = self._queue(citation_id)
        self._fetch_queue()
        return entry

    def _load_citation_stubs(self):
        """
        Provides Citation objects to the FetchEntries but doesn't populate them (some might even be new Citations
        not yet saved to the database).
        """
        if fetch_needing_citations := [fetch for fetch in self.id_to_fetch.values() if not fetch.citation]:
            # assume we request citations much more than we create them
            citations = Citation.objects.filter(
                id__in=[fetch.normalised_id.full_id for fetch in fetch_needing_citations]
            )
            citations_by_id: Dict[str, Citation] = {citation.pk: citation for citation in citations}
            for fetch in fetch_needing_citations:
                if existing := citations_by_id.get(fetch.normalised_id.full_id):
                    fetch.citation = existing
                    # Special case of having JSON populated but not the appropriate fields
                    # This will be the case if teh citation was migrated
                    if existing.data_json and not existing.last_loaded:
                        if existing.source in {CitationSource.PUBMED, CitationSource.PUBMED_CENTRAL}:
                            CitationFetchRequest._populate_from_entrez(existing, existing.data_json)
                        elif existing.source == CitationSource.NCBI_BOOKSHELF:
                            CitationFetchRequest._populate_from_nbk(existing, existing.data_json)
                        existing.last_loaded = now()
                        existing.save()
                        fetch.fetched = True
                else:
                    # doesn't save the citation yet, will do that when we retrieve the data
                    fetch.citation = fetch.normalised_id.get_or_create()

    @property
    def _all_citation_fetches(self):
        return list(self.id_to_fetch.values())

    @property
    def _all_citations_needing_loading(self):
        return [fetch for fetch in self._all_citation_fetches if fetch.should_refresh(self.refresh_before)]

    def _fetch_queue(self):
        self._load_citation_stubs()

        if ids_require_fetching := [fetch.normalised_id for fetch in self._all_citations_needing_loading]:
            citations_by_source = sorted(ids_require_fetching, key=lambda c: c.source)
            for source, citations_ids_by_source in itertools.groupby(citations_by_source, key=lambda c: c.source):
                citations_ids_by_source = list(citations_ids_by_source)
                if source == CitationSource.PUBMED:
                    self._load_from_entrez(entrez_db=EntrezDbType.PUBMED, ids=citations_ids_by_source)
                elif source == CitationSource.PUBMED_CENTRAL:
                    self._load_from_entrez(entrez_db=EntrezDbType.PUBMED_CENTRAL, ids=citations_ids_by_source)
                elif source == CitationSource.NCBI_BOOKSHELF:
                    self._load_from_nbk(ids=citations_ids_by_source)

        # Often the error is we asked NCBI Bookshelf for a list of IDs, but it only returned data for some of the IDs
        # so if nothing has returned data for a Citation, mark it as in error
        self._mark_error_if_not_fetched(ids_require_fetching, "No record was returned for this ID")

        for fetch in self.id_to_fetch.values():
            if fetch.fetched:
                # TODO, could look into making this a bulk update
                fetch.citation.save()

    def _mark_error_if_not_fetched(self, ids: Iterable[CitationRequest], error_message: str):
        for fetch_id in ids:
            try:
                if fetch := self.id_to_fetch.get(CitationIdNormalized.normalize_id(fetch_id)):
                    if not fetch.fetched:
                        fetch.fetched = True
                        fetch.citation.last_loaded = now()
                        fetch.citation.error = error_message
            except ValueError:
                pass

    def _fetch_to_populate(self, cit_id: CitationIdNormalized) -> Optional[Citation]:
        if entry := self.id_to_fetch.get(cit_id):
            entry.fetched = True
            entry.citation.blank_out()
            entry.citation.last_loaded = now()
            return entry.citation
        return None

    def _load_from_entrez(self, entrez_db: EntrezDbType, ids: List[CitationIdNormalized]):
        """
        Requests citation data from Entrez, populates via fetch _fetch_to_populate
        :param entrez_db: Must be EntrezDb.PMID or EntrezDb.PubMed
        :param ids: A list of IDs, all sources should match the equivilant EntrezDbType
        """
        request_ids = [citation_id.index for citation_id in ids]
        try:
            handle = Entrez.efetch(db=entrez_db, id=request_ids, rettype='medline', retmode='text')
            records = list(Medline.parse(handle))
            for record in records:
                normal_id: CitationIdNormalized
                if entrez_db == EntrezDbType.PUBMED:
                    normal_id = CitationIdNormalized(source=CitationSource.PUBMED, index=record.get("PMID"))
                elif entrez_db == EntrezDbType.PUBMED_CENTRAL:
                    normal_id = CitationIdNormalized(source=CitationSource.PUBMED_CENTRAL, index=record.get("PMC"))
                else:
                    raise ValueError(f"Unsupported EntrezDbType {entrez_db}")

                if citation := self._fetch_to_populate(normal_id):
                    CitationFetchRequest._populate_from_entrez(citation, record)

        except Exception as ex:
            # if this fails it's probably because a single id in ids ruined it for everybody
            report_exc_info(f'Error when attempting to Entrez.efetch ids {ids}')
            self._mark_error_if_not_fetched(ids, f'Error when attempting to Entrez.efetch ids {request_ids} : {str(ex)}')

    @staticmethod
    def _populate_from_entrez(citation: Citation, record: JsonObjType):
        citation.data_json = record
        citation.title = record.get("TI")
        citation.journal = record.get("SO")
        citation.journal_short = record.get("TA", citation.journal)

        # TODO could we just store published date and extract year?
        citation.year = get_year_from_date(record.get("DP"))
        if authors_list := record.get("FAU"):
            first_author = authors_list[0]
            first_author_last = first_author.split(",")[0]
            citation.authors_short = first_author_last
            citation.authors = ", ".join(authors_list)
        citation.abstract = record.get("AB")

    def _load_from_nbk(self, ids: Iterable[CitationIdNormalized]):
        """
        Requests citation data from Entrez, populates via fetch _fetch_to_populate
        :param ids: A list of IDs, all sources should be BookShelfID
        """
        for bookshelf_rid in ids:
            try:
                handle = Entrez.esearch(db=EntrezDbType.BOOKSHELF, term=bookshelf_rid.index, retmax=50)
                search_results = Entrez.read(handle)
                if id_list := search_results['IdList']:
                    handle = Entrez.esummary(db=EntrezDbType.BOOKSHELF, id=','.join(id_list))
                    results = Entrez.read(handle)

                    for record in results:
                        rid = record["RID"]
                        if CitationFetchRequest.TOP_LEVEL_NBK_RE.fullmatch(rid):
                            # result from esummary includes the book and chapters, with no real visible rhyme or reason
                            if citation := self._fetch_to_populate(CitationIdNormalized.from_parts(source=CitationSource.NCBI_BOOKSHELF, index=rid)):
                                CitationFetchRequest._populate_from_nbk(citation, record)

            except RuntimeError as run_error:
                report_exc_info(extra_data={'bookshelf_rid': bookshelf_rid})
                self._mark_error_if_not_fetched([bookshelf_rid], f'Error when attempting to Entrez.efetch ids {bookshelf_rid} : {str(run_error)}')

    @staticmethod
    def _populate_from_nbk(citation: Citation, record: JsonObjType):
        publish_date = record.get("DP") or record.get("PubDate")
        year = get_year_from_date(publish_date)
        book = record.get("Book")
        journal = f"Book Title: {book}"
        journal_short = book

        citation.data_json = record
        citation.title = record.get("Title")
        citation.journal = journal
        citation.journal_short = journal_short
        citation.year = year
