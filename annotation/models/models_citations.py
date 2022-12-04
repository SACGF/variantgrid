import datetime
import itertools
from dataclasses import dataclass, field
from datetime import timedelta
from enum import Enum
from typing import Optional, Iterable, List, Dict, Set, Union

from Bio import Entrez, Medline
from django.db import models
from django.db.models import TextField
from django.utils.timezone import now
from django_extensions.db.models import TimeStampedModel

from annotation.models import CitationSource
from library.log_utils import report_message, report_exc_info
from library.utils import JsonObjType
import re


def get_year_from_date(date_published) -> str:
    year = None
    if date_published:
        year = date_published.split()[0]
    return year


class CitationSource2(models.TextChoices):
    PUBMED = 'PMID', 'PMID'
    NCBI_BOOKSHELF = 'NBK', 'NCBIBookShelf'
    PUBMED_CENTRAL = 'PubMedCentral', 'PubMedCentral'


class EntrezDbType(str, Enum):
    PUBMED = "pubmed"
    PUBMED_CENTRAL = "pmc"


class Citation2(TimeStampedModel):
    # details about the ID of the citation
    id = TextField(primary_key=True)

    citation_source = models.CharField(max_length=15, choices=CitationSource.choices)
    index = models.TextField(null=True)

    # details about retrieving the JSON from a source
    data_json = models.JSONField(null=True, blank=True)
    last_loaded = models.DateTimeField(null=True, blank=True)
    error = models.TextField(null=True, blank=True)

    # details about the citation directly
    title = models.TextField(blank=True)
    journal = models.TextField(blank=True)
    journal_short = models.TextField(blank=True)
    year = models.TextField(blank=True)
    authors = models.TextField(blank=True)
    authors_short = models.TextField(blank=True)
    abstract = models.TextField(blank=True)

    def blank_out(self):
        self.title = ''
        self.journal = ''
        self.journal_short = ''
        self.year = ''
        self.authors = ''
        self.authors_short = ''
        self.abstract = ''


@dataclass(frozen=True)
class CitationIdNormalized:
    source: CitationSource2
    index: str

    @property
    def full_id(self):
        if self.source == CitationSource2.PUBMED:
            return f"PMID:{self.index}"
        elif self.source == CitationSource2.PUBMED_CENTRAL:
            return f"PMCID:{self.index}"
        else:
            return f"NBK{self.index}"

    CITATION_SPLIT_RE = re.compile(r"(?P<prefix>[a-z])+\s*(?P<semicolon>:)?\s*(?P<number_prefix>[a-z]+)?(?P<number>[0-9]+)", re.IGNORECASE)

    @staticmethod
    def from_citation(citation: Citation2):
        return CitationIdNormalized(source=citation.source, index=citation.index)

    @staticmethod
    def normalize_id(citation_id: Union[str, 'CitationIdNormalized']) -> 'CitationIdNormalized':
        if isinstance(citation_id, CitationIdNormalized):
            return citation_id
        citation_id = citation_id.strip().upper()
        if parts := Citation2.CITATION_SPLIT_RE.match():
            prefix = parts.group('prefix')
            number = parts.group('number')
            if prefix == "PUBMED":
                prefix = "PMID"
            elif prefix == "PUBMEDCENTRAL":
                prefix = "PubMedCentral"
            source = CitationSource2(prefix)
            if source == CitationSource2.PUBMED_CENTRAL:
                number = f"PMC{number}"
            return CitationIdNormalized(source=source, index=number)
        raise ValueError(f"Unable to parse Citation ID {citation_id}")


@dataclass
class CitationFetchEntry:
    normalised_id: Optional[CitationIdNormalized] = None
    citation: Optional[Citation2] = None
    requested_ids: Set[str] = field(default_factory=set)
    fetched: bool = False

    def as_citation(self):
        return Citation2(
            id=self.normalised_id,
            source=self.normalised_id.source,
            index=self.normalised_id.index
        )

    def should_refersh(self, refresh_before: datetime) -> bool:
        if last_loaded := self.citation.last_loaded:
            return last_loaded < refresh_before
        else:
            return True

    @property
    def is_valid(self) -> bool:
        return self.citation is not None

    @staticmethod
    def requested_ids_to_fetch(citation_ids: Iterable[str]) -> List['CitationFetchEntry']:
        normal_to_fetch: Dict[CitationIdNormalized, CitationFetchEntry] = dict()
        invalid_ids: Set[str] = set()

        for citation_id in citation_ids:
            try:
                normalized = CitationIdNormalized(citation_id)
                existing = normal_to_fetch.get(normalized)
                if not existing:
                    entry = CitationFetchEntry(
                        normalised_id=normalized,
                        citation=Citation2.objects.get_or_create(
                            id=normalized.full_id,
                            citation_source=normalized.source,
                            index=normalized.index
                        )
                    )
                    normal_to_fetch[normalized] = entry
                entry.requested_ids.add(citation_id)
            except:
                invalid_ids.append(citation_id)

        not_created_yet = set(normal_to_fetch.keys())
        fetch_these = set(not_created_yet)

        # assume that for any citation, it's not going to exist once, then it is going to exist 20 times
        # so try fetching first
        for citation in Citation2.objects.filter(id__in=[f.normalized_id.full_id for f in fetch_these]):
            fetched_id = CitationIdNormalized.from_citation(citation)
            not_created_yet.remove(fetched_id)
            normal_to_fetch[fetched_id].citation = citation

        if not_created_yet:
            # there were some citations that aren't already in the database
            # create them, then populate the fetch requests as if they were always tehre
            Citation2.objects.bulk_create(
                objs=[not_yet.as_citation() for not_yet in not_created_yet]
            )
            for citation in Citation2.objects.filter(id__in=[f.normalized_id.full_id for f in not_created_yet]):
                fetched_id = CitationIdNormalized.from_citation(citation)
                not_created_yet.remove(fetched_id)
                normal_to_fetch[fetched_id].citation = citation

        return normal_to_fetch.values() + [CitationFetchEntry(requested_ids=set([citation_id]))]


class CitationFetchRequest:

    def __init__(self, fetch_citations: List[CitationFetchEntry], max_age: Optional[timedelta] = None):
        requires_reloading_before = now() - (max_age if max_age else timedelta(years=100))
        requires_loading = [fetch for fetch in fetch_citations if fetch.should_refersh(requires_reloading_before)]

        self.id_to_fetch: Dict[CitationIdNormalized, CitationFetchEntry] = {fc.normalised_id: fc for fc in requires_loading}

        if requires_loading:
            citations_by_source = sorted([fetch.normalised_id for fetch in fetch_citations], lambda c: c.ontology_service)
            for source, citations_ids_by_source in itertools.groupby(citations_by_source, key=lambda c: c.ontology_service):
                citations_ids_by_source = list(citations_ids_by_source)
                if source == CitationSource2.PUBMED:
                    self.load_from_entrez(entrez_db=EntrezDbType.PUBMED, ids=citations_ids_by_source)
                elif source == CitationSource2.PUBMED_CENTRAL:
                    self.load_from_entrez(entrez_db=EntrezDbType.PUBMED_CENTRAL, ids=citations_ids_by_source)
                elif source == CitationSource2.NCBI_BOOKSHELF:
                    self.load_from_nbk(ids=citations_ids_by_source)

        for fetch in self.id_to_fetch.values():
            # would like to bulk update but that wouldn't update modified date
            fetch.citation.save()

    def mark_error_if_not_fetched(self, ids: Iterable[str, CitationIdNormalized], error_message: str):
        for fetch_id in ids:
            if fetch := self.id_to_fetch.get(fetch_id):
                fetch.loaded = True
                fetch.error = error_message

    def _fetch_to_populate(self, id: Union[str, CitationIdNormalized]) -> Optional[Citation2]:
        if entry := self.id_to_fetch.get(CitationIdNormalized.normalize_id(id)):
            entry.fetched = True
            entry.citation.blank_out()
            entry.citation.last_loaded = now()
            return entry.citation
        return None

    def load_from_entrez(self, entrez_db: EntrezDbType, ids: List[CitationIdNormalized]):
        request_ids = [id.index for id in ids]
        try:
            handle = Entrez.efetch(db=entrez_db, id=request_ids, rettype='medline', retmode='text')
            records = list(Medline.parse(handle))
            for record in records:
                normal_id: CitationIdNormalized
                if entrez_db == EntrezDbType.PUBMED:
                    normal_id = CitationIdNormalized(source=CitationSource2.PUBMED, index=record.get("PMID"))
                elif entrez_db == EntrezDbType.PUBMED_CENTRAL:
                    normal_id = CitationIdNormalized(source=CitationSource2.PUBMED_CENTRAL, index=record.get("PMC"))

                if citation := self._fetch_to_populate(normal_id):
                    CitationFetchRequest.populate_from_entrez(citation, record)

        except Exception as ex:
            # if this fails it's probably because a single id in ids ruined it for everybody
            report_exc_info(f'Error when attempting to Entrez.efetch ids {ids}')
            self.mark_error_if_not_fetched(ids, f'Error when attempting to Entrez.efetch ids {request_ids} : {str(ex)}')

    @staticmethod
    def populate_from_entrez(citation: Citation2, record: JsonObjType):
        citation.title = record.get("TI")
        citation.journal = record.get("SO")
        citation.journal_short = record.get("TA", citation.journal)
        # TODO could year be an int?
        # or could we just store published date and extract year?
        citation.year = get_year_from_date(record.get("DP"))
        if authors_list := record.get("FAU"):
            first_author = authors_list[0]
            first_author_last = first_author.split(",")[0]
            citation.authors_short = first_author_last
            citation.authors = ", ".join(authors_list)
        citation.abstract = record.get("AB")

    def load_from_nbk(self, ids: Iterable[CitationIdNormalized]):
        for bookshelf_rid in ids:
            try:
                handle = Entrez.esearch(db="books", term=bookshelf_rid.full_id, retmax=20)
                search_results = Entrez.read(handle)
                handle = Entrez.esummary(db="books", id=','.join(search_results['IdList']))
                results = Entrez.read(handle)

                for record in results:
                    if citation := self._fetch_to_populate(record["RID"]):
                        CitationFetchRequest.populate_from_nbk(citation, record)

            except RuntimeError as re:
                report_message('Searching for bookshelf_rid caused an error', level='error',
                               extra_data={'bookshelf_rid': bookshelf_rid, 'error': str(re)})
                self.mark_error_if_not_fetched([bookshelf_rid],
                                               f'Error when attempting to Entrez.efetch ids {bookshelf_rid} : {str(re)}')

    @staticmethod
    def populate_from_nbk(citation: Citation2, record: JsonObjType):
        year = get_year_from_date(record.get("DP"))
        book = record.get("Book")
        journal = f"Book Title: {book}"
        journal_short = book

        citation.title = record.get("Title")
        citation.journal = journal
        citation.journal_short = journal_short
        citation.year = year
