from dataclasses import dataclass, field
from datetime import timedelta, datetime
from functools import cached_property
from urllib.error import HTTPError
import time
from django.utils import timezone
from typing import Optional, List
import re
import json
from Bio import Entrez
from django.db import transaction

from annotation.models import ClinVarRecord, ClinVarRecordCollection
from library.log_utils import report_message
from library.utils.xml_utils import XmlParser, parser_path, PP

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
class ClinVarFetchResponse:
    """
    This warps a clinvar_record_collection just so we can override records to match the number of stars
    requested
    """

    clinvar_record_collection: ClinVarRecordCollection
    min_stars: int

    @property
    def clinvar_variation_id(self):
        return self.clinvar_record_collection.clinvar_variation_id

    @property
    def rcvs(self):
        return self.clinvar_record_collection.rcvs

    @property
    def last_loaded(self):
        return self.clinvar_record_collection.last_loaded

    @cached_property
    def records(self):
        # it's possible that we've retrieved records with min stars of 1
        # and now we re-used the cache, but only want records with min stars of 3 (so filter the others out)
        return self.clinvar_record_collection.records_with_min_stars(self.min_stars)


@dataclass
class ClinVarFetchRequest:
    """
    Is used to retrieve individual ClinVar records from the ClinVar web service for a given clinvar_variation_id
    """

    clinvar_variation_id: int
    min_stars: int
    """
    The API currently only allows us to ask for all records for that clinvar_variation_id and then we can cache
    only the ones we want. if we've already cached records (but with a lower min_stars, we can just re-use that)
    """

    max_cache_age: timedelta = field(default=timedelta(days=CLINVAR_RECORD_CACHE_DAYS))
    """
    How old until the cache is considered stale, provide seconds=0 if you want to force a refresh
    """

    def fetch(self) -> ClinVarFetchResponse:
        fetch_date = timezone.now()

        with transaction.atomic():
            clinvar_record_collection, created = ClinVarRecordCollection.objects.get_or_create(clinvar_variation_id=self.clinvar_variation_id)
            # the select_for_update() stops two simultaneous requests for updating the same clinvar_record_collection
            clinvar_record_collection = ClinVarRecordCollection.objects.select_for_update().get(pk=clinvar_record_collection.pk)
            wipe_old_records = False
            fetch_from_clinvar = True
            if not created:
                if \
                        (clinvar_record_collection.parser_version == ClinVarXmlParser.PARSER_VERSION) and \
                        clinvar_record_collection.min_stars_loaded is not None and \
                        (clinvar_record_collection.min_stars_loaded <= self.min_stars) and \
                        clinvar_record_collection.last_loaded and \
                        ((fetch_date - clinvar_record_collection.last_loaded) <= self.max_cache_age):
                    # if all the above is true, then our cache is fine
                    fetch_from_clinvar = False
                else:
                    wipe_old_records = True

            if fetch_from_clinvar:
                # so while Entrez does automatically retry on 500s, ClinVar has been providing 400s (Bad Request) when
                # the request is fine
                attempt_count = 2
                while True:
                    # loop is broken out of if it works, or raise if it fails after attempt_count
                    try:
                        response = ClinVarXmlParser.load_from_clinvar_id(clinvar_record_collection)

                        # update our cache
                        clinvar_record_collection.last_loaded = fetch_date
                        clinvar_record_collection.min_stars_loaded = self.min_stars
                        clinvar_record_collection.rcvs = response.rcvs
                        clinvar_record_collection.parser_version = ClinVarXmlParser.PARSER_VERSION
                        clinvar_record_collection.save()

                        if wipe_old_records:
                            # We *could* try to update based on SCV, and delete missing records / insert other records
                            # but a wipe and replace is easier
                            ClinVarRecord.objects.filter(clinvar_record_collection=clinvar_record_collection).delete()

                        min_star_records = [r for r in response.all_records if r.stars >= self.min_stars]
                        ClinVarRecord.objects.bulk_create(min_star_records)
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

        return ClinVarFetchResponse(
            clinvar_record_collection=clinvar_record_collection,
            min_stars=self.min_stars
        )


@dataclass
class ClinVarXmlParserOutput:
    """
    When we retrieve from ClinVar, this is our result.
    """
    rcvs: List[str]
    """
    rcvs Are condition/variant combinations - a single variant might have multiple of these to retrieve.
    The list of rcvs is only useful for debugging purposes.
    """
    all_records: List[ClinVarRecord]
    """
    All records that ClinVar had for the a given clinvar_variation_id regardless of stars.
    Up to the caller to only save the relevant records.
    Note that these will always be freshly created ClinVarRecords (which may or may not clash with existing records
    unless the caller wipes data first).
    """


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
    def load_from_clinvar_id(clinvar_record_collection: ClinVarRecordCollection) -> ClinVarXmlParserOutput:
        """
        :param clinvar_record_collection: The ClinVarRecordCollection the records should link to, also provides
        the clinvar_variation_id for us to query on.
        """

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
                            # rcv is a combination of a condition and variant of which there's 1 or more
                            # we need to retrieve all of them to get the full data
                            all_rcvs = rcvs

        if all_rcvs:
            handle = Entrez.efetch(db="clinvar", rettype="clinvarset", id=all_rcvs)
            parsed_results = ClinVarXmlParser.load_from_input(handle, clinvar_record_collection=clinvar_record_collection)
            handle.close()

        return ClinVarXmlParserOutput(
            rcvs=all_rcvs,
            all_records=parsed_results
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

            if final_value:
                self.latest.condition = final_value
