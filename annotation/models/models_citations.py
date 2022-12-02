import datetime
import itertools
from dataclasses import dataclass, field
from datetime import timedelta
from enum import Enum
from typing import Optional, Iterable, List, Dict, Set

from Bio import Entrez, Medline
from django.db import models
from django.db.models import TextField
from django.utils.timezone import now
from django_extensions.db.models import TimeStampedModel

from annotation.models import CitationSource
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
    citation_id = models.TextField(blank=True)
    citation_link = models.TextField(blank=True)
    source = models.TextField(blank=True)
    abstract = models.TextField(blank=True)

    @staticmethod
    def ensure_populated(citations: Iterable['Citation2'], max_age: Optional[timedelta] = None):
        bookshelf: List[Citation2] = []
        pmids: List[Citation2] = []

        requires_reloading_before = now() - (max_age if max_age else timedelta(years=100))
        require_loading = [cit for cit in citations if not cit.last_loaded or cit.last_loaded <= requires_reloading_before]

        citations_by_source = sorted(require_loading, lambda c: c.ontology_service)
        for source, citations_by_source in itertools.groupby(citations_by_source, key=lambda c: c.ontology_service):
            if source == CitationSource2.NCBI_BOOKSHELF:
                bookshelf += list(citations_by_source)
            else:
                pmids += list(citations_by_source)

        # now load from external source

    def update_from_data_json(self):
        self.title = ''
        self.journal = ''
        self.journal_short = ''
        self.year = ''
        self.authors = ''
        self.authors_short = ''
        self.citation_id = ''
        self.citation_link = ''
        self.source = ''
        self.abstract = ''

        if record := self.data_json:
            # PMID, PubMed Central record
            if citation_id := record.get("PMID"):
                self.title = record.get("TI")
                self.journal = record.get("SO")
                self.journal_short = record.get("TA", self.journal)
                # TODO could year be an int?
                # or could we just store published date and extract year?
                self.year = get_year_from_date(record.get("DP"))
                if authors_list := record.get("FAU"):
                    first_author = authors_list[0]
                    first_author_last = first_author.split(",")[0]
                    self.authors_short = first_author_last
                    self.authors = ", ".join(authors_list)
                self.citation_id = citation_id
                self.citation_link = f"https://www.ncbi.nlm.nih.gov/pubmed/{citation_id}"
                self.source=self.get_ontology_source_display()
                self.abstract=record.get("AB")

                # NCBI Bookshelf
            elif citation_id := record.get("RID"):
                book = record.get("Book")
                self.title = record.get("Title")
                self.journal = f"Book Title: {book}"
                self.journal_short = book
                self.year = get_year_from_date(record.get("DP"))
                self.citation_id = citation_id
                self.citation_link = f"https://www.ncbi.nlm.nih.gov/books/{citation_id}"
                self.source="NCBIBookShelf"
                #authors, authors_short, abstract = None


@dataclass(frozen=True)
class CitationIdNormalized:
    source: CitationSource2
    index: str

    @property
    def full_id(self):
        if self.source == CitationSource2.PUBMED:
            return f"PMID:{self.index}"
        elif self.source == CitationSource2.PUBMED_CENTRAL:
            raise ValueError("FIXME need to know proper layout for PubMedCentral")
        else:
            return f"NBK{self.index}"

    CITATION_SPLIT_RE = re.compile("(?P<prefix>[a-z])+(?P<semicolon>:)?(?P<number>[0-9]+)", re.IGNORECASE)

    @staticmethod
    def from_citation(citation: Citation2):
        return CitationIdNormalized(source=citation.source, index=citation.index)

    @staticmethod
    def normalize_id(citation_id: str) -> 'CitationIdNormalized':
        citation_id = citation_id.strip().upper()
        if parts := Citation2.CITATION_SPLIT_RE.match():
            prefix = parts.group('prefix')
            number = parts.group('number')
            if prefix == "PUBMED":
                prefix = "PMID"
            elif prefix == "PUBMEDCENTRAL":
                prefix = "PubMedCentral"
            source = CitationSource2(prefix)
            return CitationIdNormalized(source=source, index=number)
        raise ValueError(f"Unable to parse Citation ID {citation_id}")


@dataclass
class CitationFetchEntry:
    normalised_id: Optional[CitationIdNormalized] = None
    citation: Optional[Citation2] = None
    requested_ids: Set[str] = field(default_factory=set)

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
        if requires_loading:
            citations_by_source = sorted([fetch.citation for fetch in fetch_citations], lambda c: c.ontology_service)
            by_source: Dict[CitationSource2, List[Citation2]] = dict()
            for source, citations_by_source in itertools.groupby(citations_by_source, key=lambda c: c.ontology_service):
                by_source[source] = list(citations_by_source)




    @staticmethod
    def create_cached_citations_from_entrez(entrez_db: str, ids: Iterable[str]) -> Dict[str, JsonObjType]:
        try:
            handle = Entrez.efetch(db=entrez_db, id=ids, rettype='medline', retmode='text')
            records = list(Medline.parse(handle))
            if records:


                for cvc, record in zip(cvcs_to_query, records):
                    cc = cache_citation(cvc, record)
                    try:
                        citations_by_cvc_id[cvc.pk] = get_citation_from_cached_citation(cc)
                    except CitationException:
                        citations_by_cvc_id[cvc.pk] = _citation_error(cc)
            else:
                for cvc in cvcs_to_query:
                    citations_by_cvc_id[cvc.pk] = _citation_error(cvc)

        except Exception:
            # if this fails it's probably because a single id in ids ruined it for everybody
            for cvc in cvcs_to_query:
                citations_by_cvc_id[cvc.pk] = _citation_error(cvc)
                report_exc_info(extra_data={f'Error when attempting to Entrez.efetch ids {ids}'})

        return citations_by_cvc_id