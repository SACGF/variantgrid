from dataclasses import dataclass, field
from datetime import timedelta, datetime
from urllib.error import HTTPError
import time
from django.utils import timezone
from enum import Enum
from typing import Optional, List
import re
import json
from Bio import Entrez
from django.db import transaction

from annotation.models import ClinVarRecord, ClinVarRecordCollection
from annotation.utils.clinvar_constants import CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE
from library.log_utils import report_message
from library.utils.xml_utils import XmlParser, parser_path, PP

CLINVAR_RECORD_CACHE_DAYS = 60

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
    # FIXME, need to confirm this is 4 stars
    "practise guidelines": 4
}


@dataclass
class ClinVarFetchRequest:
    clinvar_variation_id: int
    min_stars: int
    max_cache_age: timedelta = field(default=timedelta(days=CLINVAR_RECORD_CACHE_DAYS))

    def fetch(self) -> ClinVarRecordCollection:

        fetch_date = timezone.now()

        with transaction.atomic():
            clinvar_record_collection, created = ClinVarRecordCollection.objects.get_or_create(clinvar_variation_id=self.clinvar_variation_id)
            # stops two simultaneous requests for updating the same clinvar_record_collection
            clinvar_record_collection = ClinVarRecordCollection.objects.select_for_update().get(pk=clinvar_record_collection.pk)
            wipe_old_records = False
            if not created:
                if \
                        (clinvar_record_collection.parser_version == ClinVarXmlParser.PARSER_VERSION) and \
                        clinvar_record_collection.min_stars_loaded is not None and \
                        (clinvar_record_collection.min_stars_loaded <= self.min_stars) and \
                        clinvar_record_collection.last_loaded and \
                        ((fetch_date - clinvar_record_collection.last_loaded) <= self.max_cache_age):
                    # if all the above is true, then our cache is fine

                    return clinvar_record_collection
                else:
                    wipe_old_records = True

            attempt_count = 2
            while True:
                try:
                    response = ClinVarXmlParser.load_from_clinvar_id(clinvar_record_collection).with_min_stars(self.min_stars)

                    # update our cache
                    clinvar_record_collection.last_loaded = fetch_date
                    clinvar_record_collection.min_stars_loaded = self.min_stars
                    clinvar_record_collection.rcvs = response.rcvs
                    clinvar_record_collection.parser_version = ClinVarXmlParser.PARSER_VERSION
                    clinvar_record_collection.save()

                    if wipe_old_records:
                        ClinVarRecord.objects.filter(clinvar_record_collection=clinvar_record_collection).delete()

                    ClinVarRecord.objects.bulk_create(response.records)
                    break

                except HTTPError as http_ex:
                    if http_ex.code == 400:
                        attempt_count -= 1
                        if attempt_count > 0:
                            report_message(f"400 from Entrez when fetching ClinVarRecord, waiting then trying again", level='warning')
                            time.sleep(10)
                            continue
                    # out of attempts or not 400
                    raise

        return clinvar_record_collection


@dataclass
class ClinVarFetchResponse:
    rcvs: List[str]
    records: List[ClinVarRecord]

    def with_min_stars(self, min_stars: int) -> 'ClinVarFetchResponse':
        return ClinVarFetchResponse(
            rcvs=self.rcvs,
            records=[r for r in self.records if r.stars >= min_stars]
        )


class ClinVarRetrieveMode(str, Enum):
    EXPERT_PANEL_ONLY = "expect"
    ALL_RECORDS = "all"

    @property
    def min_stars(self):
        if self == ClinVarRetrieveMode.EXPERT_PANEL_ONLY:
            return CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE
        else:
            return 0


class ClinVarXmlParser(XmlParser):

    RE_DATE_EXTRACTOR = re.compile("([0-9]+-[0-9]+-[0-9]+).*")
    RE_GOOD_CHGVS = re.compile("^(N._[0-9]+[.][0-9]+:c[.][0-9_a-zA-Z>]+)( .*)?$")
    RE_ORPHA = re.compile("ORPHA([0-9]+)")
    PARSER_VERSION = 1  # if we start caching these in the database, this is useful to know

    @staticmethod
    def parse_xml_date(text: str) -> Optional[datetime]:
        if match := ClinVarXmlParser.RE_DATE_EXTRACTOR.match(text):
            relevant_text = match.group(1)
            return datetime.strptime(relevant_text, "%Y-%m-%d")
        return None

    @staticmethod
    def load_from_clinvar_id(clinvar_record_collection: ClinVarRecordCollection) -> ClinVarFetchResponse:

        cv_handle = Entrez.esummary(db="clinvar", retmode="json", id=clinvar_record_collection.clinvar_variation_id)
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
            parsed_results = ClinVarXmlParser.load_from_input(handle, clinvar_record_collection=clinvar_record_collection)
            handle.close()

        return ClinVarFetchResponse(
            rcvs=all_rcvs,
            records=parsed_results
        )

    @staticmethod
    def load_from_input(handle, clinvar_record_collection: ClinVarRecordCollection) -> List[ClinVarRecord]:
        parsed_results = []
        for result in ClinVarXmlParser(clinvar_record_collection=clinvar_record_collection).parse(handle):
            parsed_results.append(result)
        parsed_results.sort(reverse=True)
        return parsed_results

    def __init__(self, clinvar_record_collection: ClinVarRecordCollection):
        self.clinvar_record_collection = clinvar_record_collection
        self.latest: Optional[ClinVarRecord] = None
        super().__init__()

    def reset(self):
        if self.latest:
            self.set_yieldable(self.latest)
        self.latest = ClinVarRecord(clinvar_record_collection=self.clinvar_record_collection)

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
        if submitter_date := elem.get("submitterDate"):
            self.latest.submitter_date = ClinVarXmlParser.parse_xml_date(submitter_date)

    @parser_path(
        "ClinVarResult-Set",
        "ClinVarSet",
        "ClinVarAssertion",
        "ClinicalSignificance",
        "ReviewStatus")
    def parse_review_status(self, elem):
        review_status = elem.text
        self.latest.review_status = review_status
        self.latest.stars = CLINVAR_REVIEW_STATUS_TO_STARS.get(review_status, 0)

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
        cs = elem.text
        self.latest.clinical_significance = CLINVAR_TO_VG_CLIN_SIG.get(cs, cs)

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
        if last_evaluated := elem.get("DateLastEvaluated"):
            self.latest.date_last_evaluated = ClinVarXmlParser.parse_xml_date(last_evaluated)

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
