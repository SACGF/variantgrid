from typing import List, Optional, Dict
from annotation.citations import CitationDetails, get_citations
from annotation.models import Citation
from annotation.regexes import DbRegexes, DbRefRegexes
from library.utils import JsonObjType

"""
In future I'm looking to simplify Citation -> CachedCitation -> CitationDetails into JUST citation.
Hoping that use of CitationLoader will make things easier to transition between.
"""


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

    CITATION_SEARCH = DbRefRegexes(regexes=[DbRegexes.PUBMED, DbRegexes.NCBIBookShelf])

    def __init__(self):
        self.citation_receipts: Dict[str, CitationReceipt] = {}

    def search_ids(self, citation_str: str) -> List[CitationReceipt]:
        """
        :param citation_str: Text like "Find the citations in PMID:234343, 46634534 but don't forget NBK2343433"
        :return: A list of CitationReceipts that will be populated once load is called
        """
        receipts = list()
        for search_result in CitationLoader.CITATION_SEARCH.search(citation_str):
            citation: Citation = Citation.objects.get(pk=search_result.internal_id)
            if receipt := self.add_citation(citation, only_receipt_if_new=True):
                receipts.append(receipt)
        return receipts

    def add_dbrefs(self, db_refs: List[JsonObjType]) -> List[CitationReceipt]:
        """
        :param db_refs: An array of db_refs (complete with 'internal_id')
        :return: A list of CitationReceipts that will be populated once load is called
        """
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
        """
        Actually populates the citation receipts
        :return: All the unique CitationReceipts that have been added using the above methods
        """
        if self.citation_receipts:
            for detail in get_citations((receipt.citation for receipt in self.citation_receipts.values())):
                if receipt := self.citation_receipts.get(detail.citation_id):
                    receipt.set_loaded(detail)
            return list(self.citation_receipts.values())
        return []
