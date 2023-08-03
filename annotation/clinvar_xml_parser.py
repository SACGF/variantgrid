from datetime import datetime
from dataclasses import dataclass
from enum import Enum
from functools import cached_property
from typing import Optional, List
import re
import json
from Bio import Entrez

from annotation.utils.clinvar_constants import CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE
from library.utils.xml_utils import XmlParser, parser_path, PP

CLINVAR_TO_VG_CLIN_SIG = {
    "Benign": "B",
    "Likely benign": "LB",
    "Uncertain Significance": "VUS",
    "Likely pathogenic": "LP",
    "Pathogenic": "P",
    "risk factor": "R",
    "drug response": "D"
}

CLINVAR_REVIEW_STATUS_TO_STARS = {
    "no assertion criteria provided": 1,
    "criteria provided, single submitter": 2,
    "reviewed by expert panel": 3,
    # FIXME, need to confirm this is 4 stars
    "practise guidelines": 4
}


@dataclass
class ClinVarRecord:
    record_id: Optional[str] = None
    org_id: Optional[str] = None
    genome_build: Optional[str] = None
    review_status: Optional[str] = None
    submitter: Optional[str] = None
    submitter_date: Optional[str] = None
    c_hgvs: Optional[str] = None
    variant_coordinate: Optional[str] = None
    condition: Optional[str] = None
    clinical_significance: Optional[str] = None
    gene_symbol: Optional[str] = None
    date_last_evaluated: Optional[str] = None
    interpretation_summary: Optional[str] = None
    assertion_method: Optional[str] = None
    allele_origin: Optional[str] = None
    parser_version: int = 1
    parsed_date: Optional[datetime] = None

    @property
    def _sort_order(self):
        return self.stars, self.submitter_datetime

    def __lt__(self, other):
        return self._sort_order < other._sort_order

    @property
    def submitter_datetime(self) -> datetime:
        return datetime.strptime(self.submitter_date, "%Y-%m-%d")

    def export_condition(self):
        return self.condition or "Not set"

    @cached_property
    def stars(self) -> int:
        return CLINVAR_REVIEW_STATUS_TO_STARS.get(self.review_status, 0)

    @property
    def clinical_significance_vg(self):
        return CLINVAR_TO_VG_CLIN_SIG.get(self.clinical_significance, "VUS")

    @property
    def is_good_quality(self) -> bool:
        # if not (self.c_hgvs and self.gene_symbol):
        #     return False
        # if "?" in self.c_hgvs or "[" in self.c_hgvs or " " in self.c_hgvs:
        #     return False
        return True

    @property
    def is_expert_panel_or_greater(self):
        return self.stars >= CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE


@dataclass
class ClinVarFetchResponse:
    clinvar_variation_id: int
    rcvs: List[str]
    records: List[ClinVarRecord]
    cached: bool = False


class ClinVarRetrieveMode(str, Enum):
    EXPERT_PANEL_ONLY = "expect"
    ALL_RECORDS = "all"


def fetch_clinvar_records(clinvar_variation_id, record_mode: ClinVarRetrieveMode) -> ClinVarFetchResponse:
    """
    Call to retrieve individual ClinVarRecords
    :param clinvar_variation_id:
    :param mode:
    :return:
    """
    response = ClinVarXmlParser.load_from_clinvar_id(clinvar_variation_id)
    response.records = [r for r in response.records if r.is_expert_panel_or_greater]
    return response


class ClinVarXmlParser(XmlParser):

    RE_GOOD_CHGVS = re.compile("^(N._[0-9]+[.][0-9]+:c[.][0-9_a-zA-Z>]+)( .*)?$")
    RE_ORPHA = re.compile("ORPHA([0-9]+)")
    PARSER_VERSION = 1  # if we start caching these in the database, this is useful to know

    @staticmethod
    def load_from_clinvar_id(clinvar_variation_id) -> ClinVarFetchResponse:

        cv_handle = Entrez.esummary(db="clinvar", retmode="json", id=clinvar_variation_id)
        json_data = json.loads(cv_handle.read())
        cv_handle.close()

        all_rcvs: List[str] = []
        parsed_results: List[ClinVarRecord] = []
        if result := json_data.get("result"):
            if uuids := result.get("uids"):
                if uuid_data := result.get(uuids[0]):
                    if supporting_submissions := uuid_data.get("supporting_submissions"):
                        if rcvs := supporting_submissions.get("rcv"):
                            all_rcvs = rcvs

        if all_rcvs:
            handle = Entrez.efetch(db="clinvar", rettype="clinvarset", id=all_rcvs)
            parsed_results = ClinVarXmlParser.load_from_input(handle)
            handle.close()

        return ClinVarFetchResponse(
            clinvar_variation_id=clinvar_variation_id,
            rcvs=all_rcvs,
            records=parsed_results
        )

    @staticmethod
    def load_from_input(handle) -> List[ClinVarRecord]:
        parsed_results = []
        for result in ClinVarXmlParser().parse(handle):
            parsed_results.append(result)
        parsed_results.sort(reverse=True)
        return parsed_results

    def __init__(self):
        self.latest: Optional[ClinVarRecord] = None
        super().__init__()

    def reset(self):
        if self.latest and self.latest.is_good_quality:
            self.latest.parser_version = ClinVarXmlParser.PARSER_VERSION
            self.latest.parsed_date = datetime.now()
            self.set_yieldable(self.latest)
        self.latest = ClinVarRecord()

    def finish(self):
        self.reset()

    @parser_path("ClinVarResult-Set", "ClinVarSet", "ClinVarAssertion")
    def new_record(self, elem):
        self.reset()

    @parser_path("ClinVarResult-Set", "ClinVarSet", "ClinVarAssertion", PP("ClinVarAccession", Type="SCV"))
    def record_id(self, elem):
        self.latest.record_id = elem.get("Acc")
        self.latest.org_id = elem.get("OrgID")

    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        "ObservedIn",
        "Sample",
        "Origin")
    def parse_allele_origin(self, elem):
        self.latest.allele_origin = elem.text

    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        PP("MeasureSet", Type="Variant"),
        PP("Measure", Type="Variation"),
        "AttributeSet",
        PP("Attribute", Type="HGVS"))
    def parse_c_hgvs(self, elem):
        if hgvs := elem.text:
            if match := ClinVarXmlParser.RE_GOOD_CHGVS.match(hgvs):
                hgvs = match.group(1)
                self.latest.c_hgvs = hgvs
            self.latest.c_hgvs = hgvs

    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        PP("MeasureSet", Type="Variant"),
        PP("Measure", Type="Variation"),
        "SequenceLocation")
    def parse_variant_coordinate(self, elem):
        assembly = elem.get("Assembly")
        chr = elem.get("Chr")
        start = elem.get("start")
        ref = elem.get("referenceAllele")
        alt = elem.get("alternateAllele")
        self.latest.variant_coordinate = f"{chr}:{start} {ref}>{alt} ({assembly})"


    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        "ClinVarSubmissionID")
    def parse_reviewer(self, elem):
        self.latest.submitter = elem.get("submitter")
        self.latest.submitter_date = elem.get("submitterDate")

    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        "ClinicalSignificance",
        "ReviewStatus")
    def parse_review_status(self, elem):
        self.latest.review_status = elem.text

    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        PP("MeasureSet", Type="Variant"),
        PP("Measure", Type="Variation"),
        PP("MeasureRelationship", Type="variant in gene"),
        "Symbol",
        PP("ElementValue", Type="Preferred"))
    def parse_gene_symbol(self, elem):
        self.latest.gene_symbol = elem.text

    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        PP("MeasureSet", Type="Variant"),
        PP("Measure", Type="Variation"),
        PP("MeasureRelationship", Type="variant in gene"),
        "Symbol",
        PP("ElementValue", Type="Preferred"))
    def parse_gene_symbol(self, elem):
        self.latest.gene_symbol = elem.text

    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        "ClinicalSignificance",
        "Description")
    def parse_clinical_significance(self, elem):
        self.latest.clinical_significance = elem.text

    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        "ClinicalSignificance",
        "Comment")
    def parse_interpretation_summary(self, elem):
        self.latest.interpretation_summary = elem.text

    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        "AttributeSet",
        PP("Attribute", Type="AssertionMethod"))
    def parse_assertion_method(self, elem):
        self.latest.assertion_method = elem.text

    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        "ClinicalSignificance")
    def parse_date_last_evaluated(self, elem):
        self.latest.date_last_evaluated = elem.get("DateLastEvaluated")

    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        PP("TraitSet", Type="Disease"),
        PP("Trait", Type="Disease"),
        "XRef")
    def parse_condition(self, elem):
        if not self.latest.condition:
            db = elem.get("DB")
            id = elem.get("ID")
            final_value = None
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

            if final_value:
                self.latest.condition = final_value
