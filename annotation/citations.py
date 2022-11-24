import json
from dataclasses import dataclass
from typing import Dict, Iterable, List, Any, Union

from Bio import Entrez, Medline

from annotation.models import Citation
from annotation.models.models import CachedCitation, CitationException
from annotation.models.models_enums import CitationSource
# @param id The primary key of the Citation that led to this CitationDetails
# @param citation_id The id as it's known in the external database, e.g. PMID:9003501 has a citation_id of 9003501
from library.log_utils import report_exc_info, report_message


@dataclass
class CitationDetails:
    id: int
    title: str
    journal: str
    journal_short: str
    year: str
    authors: str
    authors_short: str
    citation_id: str
    citation_link: str
    source: str
    abstract: str
    is_error: bool = False

    def __eq__(self, other):
        return self.source == other.source and self.citation_id == other.citation_id

    def __hash__(self):
        return hash(self.source) + hash(self.citation_id)

    def __lt__(self, other):
        if self.source < other.source:
            return True
        if self.source == other.source:
            return self.citation_id.rjust(10, '0') < other.citation_id.rjust(10, '0')
        return False


CITATION_COULD_NOT_LOAD_TEXT = 'Could not load citation summary'


def _citation_error(cached_citation: Union[Citation, CachedCitation]) -> CitationDetails:
    citation: Citation
    if isinstance(cached_citation, CachedCitation):
        citation = cached_citation.citation
    else:
        citation = cached_citation

    citation_source = citation.get_citation_source_display()

    return CitationDetails(id=citation.id,
                           title=f"{citation.citation_id} - {CITATION_COULD_NOT_LOAD_TEXT}",
                           journal=None,
                           journal_short=None,
                           year=None,
                           authors=None,
                           authors_short=None,
                           citation_id=citation.citation_id,
                           citation_link=None,
                           source=citation_source,
                           abstract=None,
                           is_error=True)


def get_year_from_date(date_published) -> str:
    year = None
    if date_published:
        year = date_published.split()[0]
    return year


def get_citation_from_cached_citation(cached_citation: CachedCitation) -> CitationDetails:
    try:
        record = cached_citation.get_record_or_fail()
    except CitationException:
        return _citation_error(cached_citation)

    citation_id = record["PMID"]
    authors = record.get("FAU")
    authors_short = None
    if authors:
        first_author = authors[0]
        first_author_last = first_author.split(",")[0]
        authors_short = first_author_last

        authors = ", ".join(authors)

    journal = record.get("SO")
    journal_short = record.get("TA", journal)
    year = get_year_from_date(record.get("DP"))

    return CitationDetails(id=cached_citation.citation_id,
                           title=record.get("TI"),
                           journal=journal,
                           journal_short=journal_short,
                           year=year,
                           authors=authors,
                           authors_short=authors_short,
                           citation_id=citation_id,
                           citation_link="https://www.ncbi.nlm.nih.gov/pubmed/" + citation_id,
                           source=cached_citation.citation.get_citation_source_display(),
                           abstract=record.get("AB"))


def get_ncbi_bookshelf_citation_from_cached_citation(cached_citation: CachedCitation) -> CitationDetails:
    try:
        record = cached_citation.get_record_or_fail()
    except CitationException:
        return _citation_error(cached_citation)

    BOOK_URL = "https://www.ncbi.nlm.nih.gov/books/"
    citation_id = record["RID"]
    year = get_year_from_date(record.get("DP"))
    book = record.get("Book")
    journal = f"Book Title: {book}"
    journal_short = book

    return CitationDetails(id=cached_citation.citation_id,
                           title=record.get("Title"),
                           journal=journal,
                           journal_short=journal_short,
                           year=year,
                           authors=None,
                           authors_short=None,
                           citation_id=citation_id,
                           citation_link=BOOK_URL + citation_id,
                           source='NCBIBookShelf',
                           abstract=None)


def cache_citation(citation: Citation, record: Dict[str, Any]) -> CachedCitation:
    json_string = json.dumps(record)
    has_error = 'Error occurred' in json_string
    # There may be a race condition (multiple requests going off fetching citations)
    # So use get_or_create and don't worry about it
    cc, created = CachedCitation.objects.get_or_create(citation=citation,
                                                       defaults={
                                                           'json_string': json_string,
                                                           'has_error': has_error
                                                       })
    if created:
        # optimisation so we don't have to re-parse the json_string
        cc._record = record

    return cc


def create_cached_citations_from_entrez(entrez_db, cvcs_to_query: List[Citation]) -> Dict[int, CitationDetails]:
    citations_by_cvc_id = {}
    ids = [str(cvc.citation_id) for cvc in cvcs_to_query]
    try:
        h = Entrez.efetch(db=entrez_db, id=ids, rettype='medline', retmode='text')
        records = list(Medline.parse(h))
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


def get_summary_for_bookshelf_rid(bookshelf_rid) -> Dict[str, Any]:

    handle = Entrez.esearch(db="books", term=bookshelf_rid, retmax=20)
    search_results = None
    try:
        search_results = Entrez.read(handle)
        handle = Entrez.esummary(db="books", id=','.join(search_results['IdList']))
        results = Entrez.read(handle)

        for r in results:
            if r["RID"] == bookshelf_rid:
                return r
    except RuntimeError as re:
        report_message('Searching for bookshelf_rid caused an error', level='warning', extra_data={'bookshelf_rid':bookshelf_rid, 'error': str(re)})
        return {"id": bookshelf_rid, "status": "Error occurred", "error": str(re), "reporter": "This is a variantgrid message"}

    return {"id": bookshelf_rid, "status": "Error occurred", "error": f"No RID matching {bookshelf_rid} returned", "reporter": "This is a variantgrid message"}


def create_cached_citations_from_books(cvcs_to_query_books: Iterable[Citation]) -> Dict[int, CitationDetails]:
    citations_by_cvc_id = {}
    for cvc in cvcs_to_query_books:
        record = get_summary_for_bookshelf_rid(cvc.citation_id)
        if record:
            cc = cache_citation(cvc, record)
            citations_by_cvc_id[cvc.pk] = get_ncbi_bookshelf_citation_from_cached_citation(cc)

    return citations_by_cvc_id


def get_citations(citations: Iterable[Citation]) -> List[CitationDetails]:
    citations_by_citation_id: Dict[int, CitationDetails] = {}
    cites_to_query_pubmed: List[Citation] = []
    cites_to_query_pmc: List[Citation] = []
    cites_to_query_books: List[Citation] = []

    # Get the existing ones
    for cite in citations:
        if cite.citation_source == CitationSource.NCBI_BOOKSHELF:
            try:
                citations_by_citation_id[cite.pk] = get_ncbi_bookshelf_citation_from_cached_citation(cite.cachedcitation)
            except CachedCitation.DoesNotExist:
                cites_to_query_books.append(cite)
        else:
            try:
                citations_by_citation_id[cite.pk] = get_citation_from_cached_citation(cite.cachedcitation)
            except CachedCitation.DoesNotExist:
                if cite.citation_source == CitationSource.PUBMED:
                    cites_to_query_pubmed.append(cite)
                elif cite.citation_source == CitationSource.PUBMED_CENTRAL:
                    cites_to_query_pmc.append(cite)

    if cites_to_query_pubmed:
        pubmed_cvcs = create_cached_citations_from_entrez('pubmed', cites_to_query_pubmed)
        citations_by_citation_id.update(pubmed_cvcs)

    if cites_to_query_pmc:
        pmc_cvcs = create_cached_citations_from_entrez('pmc', cites_to_query_pmc)
        citations_by_citation_id.update(pmc_cvcs)

    if cites_to_query_books:
        books_cvcs = create_cached_citations_from_books(cites_to_query_books)
        citations_by_citation_id.update(books_cvcs)

    cached_citations: List[CitationDetails] = []
    for k in sorted(citations_by_citation_id):
        cc = citations_by_citation_id[k]
        cached_citations.append(cc)

    cached_citations.sort()

    return cached_citations
