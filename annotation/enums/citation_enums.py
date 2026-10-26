"""
Citation enums: CitationSource (the kind of citation a Citation stores, with from_legacy_code mapping the old
prefixes) and EntrezDbType (which Entrez database a CitationFetchRequest queries for each source).
"""
from enum import Enum
from typing import Optional

from django.db import models


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
            "PUBMEDCENTRAL": CitationSource.PUBMED_CENTRAL,

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
