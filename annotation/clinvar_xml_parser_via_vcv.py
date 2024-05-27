import re
from typing import Optional

from Bio import Entrez

from annotation.clinvar_xml_parser import ClinVarXmlParser, ClinVarXmlParserOutput, CLINVAR_REVIEW_STATUS_TO_STARS, \
    CLINVAR_TO_VG_CLIN_SIG, SOMATIC_CLIN_SIG_VALUE
from annotation.models import ClinVarRecord
from classification.enums import AlleleOriginBucket
from library.utils.xml_utils import parser_path, PP


class ConditionList:

    def __init__(self):
        self.disease: list[str] = []
        self.findings: list[str] = []
        self.choices: list[str] = []

    @property
    def result(self) -> Optional[str]:
        if self.disease:
            return "; ".join(self.choices)
        elif self.findings:
            return "; ".join(self.choices)
        elif self.choices:
            # maybe choices should be uncertain?
            return "; ".join(self.choices)
        else:
            return None


class ClinVarXmlParserViaVCV(ClinVarXmlParser):

    PARSER_VERSION = 210  # change this whenever the parsing changes, so we know to ignore the old cache

    RE_LEGACY_CS = re.compile("^Converted during submission to (.*?).?$")

    @classmethod
    def load_from_clinvar_id(cls, clinvar_variation_id: int) -> ClinVarXmlParserOutput:
        vcv = f"{clinvar_variation_id}"
        while len(vcv) < 9:
            vcv = f"0{vcv}"
        vcv = f"VCV{vcv}"

        handle = Entrez.efetch(db="clinvar", rettype="vcv", id=vcv)
        parsed_results = cls.load_from_input(handle)
        handle.close()

        return ClinVarXmlParserOutput(
            all_records=parsed_results,
            urls=[
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&id={vcv}&rettype=vcv"
            ]
        )

    def __init__(self):
        self.allele_origin_set = set()
        self.human = False
        self.non_human = False
        self.somatic_clinical_impact = False
        self.oncogenicity_classification = False
        self.germline_classification = False
        self.condition_list = ConditionList()

        super().__init__(
            prefix=["ClinVarResult-Set", "VariationArchive", "ClassifiedRecord", "ClinicalAssertionList", "ClinicalAssertion"]
        )

    @parser_path(on_start=True)
    def new_record(self, elem):
        self.reset()
        self.allele_origin_set.clear()
        self.human = False
        self.non_human = False
        self.somatic_clinical_impact = False
        self.oncogenicity_classification = False
        self.germline_classification = False
        self.condition_list = ConditionList()

        self.latest.allele_origin_bucket = AlleleOriginBucket.UNKNOWN
        self.latest.submitter_date = ClinVarXmlParser.parse_xml_date(elem.get("SubmissionDate"))

    def reset(self):
        super().reset()

    @parser_path("ObservedInList", "ObservedIn", "Sample", "Origin")
    def allele_origin(self, elem):
        self.allele_origin_set.add(elem.text)

    @parser_path("ObservedInList", "ObservedIn", "Sample", "Species")
    def species(self, elem):
        if elem.text == "human" or ("homo" in elem.text.lower() and "sapiens" in elem.text.lower()):
            self.human = True
        else:
            self.non_human = True

    @parser_path(PP("ClinVarAccession", Type="SCV"))
    def record_id(self, elem):
        self.latest.record_id = elem.get("Accession")
        self.latest.org_id = elem.get("OrgID")
        self.latest.date_clinvar_created = ClinVarXmlParser.parse_xml_date(elem.get("DateCreated"))
        self.latest.date_clinvar_updated = ClinVarXmlParser.parse_xml_date(elem.get("DateUpdated"))
        self.latest.submitter = elem.get("SubmitterName")

    @parser_path("ClinVarSubmissionID")
    def parse_submitted_assembly(self, elem):
        self.latest.genome_build = elem.get("submittedAssembly")

    @parser_path(
        "Classification",
        "ReviewStatus")
    def parse_review_status(self, elem):
        review_status = elem.text
        self.latest.review_status = review_status
        self.latest.stars = CLINVAR_REVIEW_STATUS_TO_STARS.get(review_status, 0)

    @parser_path("Classification")
    def parse_date_last_evaluated(self, elem):
        if dle := elem.get("DateLastEvaluated"):
            self.latest.date_last_evaluated = ClinVarXmlParser.parse_xml_date(dle)

    @parser_path(
        "Classification",
        "Comment")
    def parse_interpretation_summary(self, elem):
        self.latest.interpretation_summary = elem.text

    @parser_path(
        "Classification",
        "GermlineClassification")
    def parse_clinical_significance_desc(self, elem):
        if cs := elem.text:
            self.germline_classification = True
            cs = cs.lower()
            self.latest.clinical_significance = CLINVAR_TO_VG_CLIN_SIG.get(cs, cs)

    @parser_path(
        "Classification",
        "SomaticClinicalImpact")
    def parse_somatic_clinical_significance_desc(self, elem):
        if cs := elem.text:
            self.somatic_clinical_impact = True
            cs = cs.lower()
            cs = cs.split("-")[0].strip()
            self.latest.somatic_clinical_significance = SOMATIC_CLIN_SIG_VALUE.get(cs, cs)

    @parser_path(
        "Classification",
        "OncogenicityClassification",
        "Description")
    def parse_oncogenicity_classification(self, elem):
        if cs := elem.text:
            self.oncogenicity_classification = True
            cs = cs.lower()
            self.latest.clinical_significance = CLINVAR_TO_VG_CLIN_SIG.get(cs, cs)

    @parser_path(
        "Classification",
        "OncogenicityClassification")
    def parse_oncogenicity_classification(self, elem):
        if cs := elem.text:
            self.oncogenicity_classification = True
            cs = cs.lower()
            self.latest.clinical_significance = CLINVAR_TO_VG_CLIN_SIG.get(cs, cs)

    @parser_path(
        "Classification",
        "Comment")
    def parse_clinical_significance_comment(self, elem):
        # prioritise comment if it's in the format of
        # Converted during submission to XXX
        if cs := elem.text:
            if match := ClinVarXmlParserViaVCV.RE_LEGACY_CS.match(cs):
                match_cs = match.group(1).lower()
                self.latest.clinical_significance = CLINVAR_TO_VG_CLIN_SIG.get(match_cs, match_cs)

    @parser_path(
        "SimpleAllele",
        "OtherNameList",
        "Name"
    )
    def _parse_possible_c_hgvs(self, elem):
        if not elem.get("Type"):
            self.assign_better_hgvs(elem.text)

    @parser_path(
        "SimpleAllele",
        "AttributeSet",
        PP("Attribute", Type="HGVS"))
    def _parse_c_hgvs(self, elem):
        self.assign_better_hgvs(elem.text)

    @parser_path(
        "SimpleAllele",
        "HGVSlist",
        "HGVS",
        "Expression")
    def parse_hgvs_2(self, elem):
        self.assign_better_hgvs(elem.text)

    @parser_path(
        "SimpleAllele",
        "Location",
        "SequenceLocation"
    )
    def parse_variant_coordinate(self, elem):
        assembly = elem.get("Assembly")
        chr = elem.get("Chr")
        start = elem.get("start")
        ref = elem.get("referenceAllele")
        alt = elem.get("alternateAllele")
        if chr and start and ref and alt:
            self.latest.variant_coordinate = f"{chr}:{start} {ref}>{alt} ({assembly})"

    @parser_path(
        PP("TraitSet", Type="Disease"),
        PP("Trait", Type="Disease"),
        "XRef")
    def parse_condition(self, elem):
        if not self.latest.condition or ":" not in self.latest.condition:
            if final_value := ClinVarXmlParser.clean_term_id_from_elem(elem):
                self.condition_list.disease.append(final_value)

    @parser_path(
        PP("TraitSet", Type="Finding"),
        PP("Trait", Type="Finding"),
        "XRef")
    def parse_condition_finding(self, elem):
        if not self.latest.condition or ":" not in self.latest.condition:
            if final_value := ClinVarXmlParser.clean_term_id_from_elem(elem):
                self.condition_list.findings.append(final_value)

    @parser_path(
        PP("TraitSet", Type="TraitChoice"),
        PP("Trait", Type="Disease"),
        "XRef")
    def parse_trait_choices(self, elem):
        if not self.latest.condition or ":" not in self.latest.condition:
            if final_value := ClinVarXmlParser.clean_term_id_from_elem(elem):
                self.condition_list.choices.append(final_value)

    @parser_path(
        PP("TraitSet", Type="DrugResponse"),
        PP("Trait", Type="DrugResponse"),
        "Name",
        PP("ElementValue", Type="Preferred"))
    def parse_drug_response(self, elem):
        if not self.latest.condition:
            self.condition_list.disease.append(elem.text)

    @parser_path(
        PP("TraitSet", Type="Disease"),
        PP("Trait", Type="Disease"),
        "Name",
        PP("ElementValue", Type="Preferred"))
    def parse_condition_name(self, elem):
        if not self.latest.condition:
            condition_text = elem.text
            if condition_text.lower() not in ("not specified", "not provided"):
                self.condition_list.disease.append(elem.text)

    def post_record_parse(self, obj: ClinVarRecord):
        if self.non_human and not self.human:
            obj.mark_invalid()
        else:
            allele_origins = self.allele_origin_set
            if self.somatic_clinical_impact or self.oncogenicity_classification:
                obj.allele_origin_bucket = AlleleOriginBucket.SOMATIC
            elif "germline" in allele_origins:
                obj.allele_origin_bucket = AlleleOriginBucket.GERMLINE
            elif "somatic" in allele_origins:
                obj.allele_origin_bucket = AlleleOriginBucket.SOMATIC
            else:
                obj.allele_origin_bucket = AlleleOriginBucket.UNKNOWN
        self.latest.condition = self.condition_list.result
