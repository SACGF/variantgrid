from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import datetime
from typing import Optional, List
import re

from annotation.models import ClinVarRecord
from library.utils.xml_utils import XmlParser

"""
This file is responsible for retrieving data from the ClinVar API end points to get more granular data about a given
ClinVar variant.
"""

CLINVAR_RECORD_CACHE_DAYS = 60
"""
Number of days we keep ClinVar records cached before we will re-ask ClinGen for them
"""

CLINVAR_TO_VG_CLIN_SIG = {
    "Benign": "B",
    "Likely benign": "LB",
    "Uncertain significance": "VUS",
    "Likely pathogenic": "LP",
    "Pathogenic": "P",
    "risk factor": "R",
    "drug response": "D"
}

CLINVAR_REVIEW_STATUS_TO_STARS = {
    "no assertion criteria provided": 1,
    "criteria provided, single submitter": 2,
    "reviewed by expert panel": 3,
    "practise guidelines": 4
}


@dataclass
class ClinVarXmlParserOutput:
    """
    When we retrieve from ClinVar, this is our result.
    """
    urls: List[str]
    """
    urls Are condition/variant combinations - a single variant might have multiple of these to retrieve.
    The list of urls is only useful for debugging purposes.
    """
    all_records: List[ClinVarRecord]
    """
    All records that ClinVar had for the a given clinvar_variation_id regardless of stars.
    Up to the caller to only save the relevant records.
    Note that these will always be freshly created ClinVarRecords (which may or may not clash with existing records
    unless the caller wipes data first).
    """


class ClinVarXmlParser(XmlParser, ABC):

    RE_DATE_EXTRACTOR = re.compile("([0-9]+-[0-9]+-[0-9]+).*")
    RE_GOOD_CHGVS = re.compile("^((N._|ENST)(.+)?:[acgn][.][^ ]+)( .*)?$")
    RE_ORPHA = re.compile("ORPHA([0-9]+)")

    @staticmethod
    def parse_xml_date(text: str) -> Optional[datetime]:
        if match := ClinVarXmlParser.RE_DATE_EXTRACTOR.match(text):
            relevant_text = match.group(1)
            return datetime.strptime(relevant_text, "%Y-%m-%d")
        return None

    def assign_better_hgvs(self, text: str) -> str:
        if text:
            if match := ClinVarXmlParser.RE_GOOD_CHGVS.match(text):
                text = match.group(1)
                def hgvs_score(some_hgvs):
                    if some_hgvs is None:
                        return 0
                    if not "c." in some_hgvs:
                        return 1
                    if not "(" in some_hgvs:
                        return 2
                    return 3

                if hgvs_score(text) > hgvs_score(self.latest.c_hgvs):
                    self.latest.c_hgvs = text

    @classmethod
    @abstractmethod
    def load_from_clinvar_id(cls, clinvar_variation_id: int) -> ClinVarXmlParserOutput:
        pass

    @classmethod
    def load_from_input(cls, handle) -> List[ClinVarRecord]:
        parsed_results: List[ClinVarRecord] = []
        for result in cls().parse(handle):
            parsed_results.append(result)
        parsed_results.sort(reverse=True)
        return parsed_results

    def __init__(self, prefix: Optional[List] = None):
        self.latest: Optional[ClinVarRecord] = None
        super().__init__(prefix=prefix)

    def reset(self):
        if self.latest:
            self.set_yieldable(self.latest)
        self.latest = ClinVarRecord()

    def finish(self):
        self.reset()

    @staticmethod
    def clean_term_id_from_elem(elem):
        db = elem.get("DB")
        id = elem.get("ID")
        final_value = None
        # ClinVar is always a bit weird about how it represents different conditions
        # use this to convert to standard term descriptions
        if db == "MONDO":
            final_value = id
        elif db == "OMIM":
            final_value = f"OMIM:{id}"
        elif db == "Orphanet":
            if m := ClinVarXmlParser.RE_ORPHA.match(id):
                final_value = f"ORPHA:{m.group(1)}"
        elif db == "MedGen":
            final_value = f"MedGen:{id}"
        elif db == "HP":
            final_value = id
        elif db == "GeneReviews":
            final_value = id
        elif db.upper() == "MESH":
            final_value = f"MeSH:{id}"
        else:
            final_value = f"{db} {id}"
        return final_value
