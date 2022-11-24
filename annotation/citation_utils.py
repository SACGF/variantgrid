from typing import List, Optional, Dict
from annotation.citations import CitationDetails, get_citations
from annotation.models import Citation
from annotation.regexes import DbRegexes, DbRefRegexes
from library.utils import JsonObjType


class CitationReceipt:

    def __init__(self, citation: Citation):
        self.citation = citation
        self.details: Optional[CitationDetails] = None

    def set_loaded(self, details: CitationDetails):
        self.details = details

    @property
    def is_valid(self) -> bool:
        return self.details and not self.details.is_error


class CitationLoader:
    """
    A class we can use as an intermediate if we ever simplify the
    Citation <-> CachedCitation <-> CitationDetails model
    """

    CITATION_SEARCH = DbRefRegexes(regexes=[DbRegexes.PUBMED, DbRegexes.NCBIBookShelf])

    def __init__(self):
        self.citation_receipts: Dict[str, CitationReceipt] = {}

    def search_ids(self, citation_str: str) -> List[CitationReceipt]:
        receipts = list()
        for search_result in CitationLoader.CITATION_SEARCH.search(citation_str):
            citation: Citation = Citation.objects.get(pk=search_result.internal_id)
            if receipt := self.add_citation(citation, only_receipt_if_new=True):
                receipts.append(receipt)
        return receipts

    def add_dbrefs(self, db_refs: List[JsonObjType]) -> List[CitationReceipt]:
        receipts = list()
        for db_ref in db_refs:
            if db_ref.get('db') in {DbRegexes.PUBMED.db, DbRegexes.NCBIBookShelf.db} and (internal_id := db_ref.get('internal_id')):
                citation = Citation.objects.get(pk=internal_id)
                if receipt := self.add_citation(citation, only_receipt_if_new=True):
                    receipts.append(receipt)
        return receipts

    def add_citation(self, citation: Citation, only_receipt_if_new=False) -> Optional[CitationReceipt]:
        receipt = self.citation_receipts.get(citation.citation_id)
        if not receipt:
            receipt = CitationReceipt(citation)
            self.citation_receipts[citation.citation_id] = receipt
            return receipt
        if not only_receipt_if_new:
            return receipt

    def load(self) -> List[CitationReceipt]:
        if self.citation_receipts:
            for detail in get_citations((receipt.citation for receipt in self.citation_receipts.values())):
                if receipt := self.citation_receipts.get(detail.citation_id):
                    receipt.set_loaded(detail)
            return list(self.citation_receipts.values())
        return []
