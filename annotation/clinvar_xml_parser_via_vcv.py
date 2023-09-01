from Bio import Entrez

from annotation.clinvar_xml_parser import ClinVarXmlParser, ClinVarXmlParserOutput, CLINVAR_REVIEW_STATUS_TO_STARS, \
    CLINVAR_TO_VG_CLIN_SIG
from library.utils.xml_utils import parser_path, PP


class ClinVarXmlParserViaVCV(ClinVarXmlParser):

    PARSER_VERSION = 200  # change this whenever the parsing changes, so we know to ignore the old cache

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
        super().__init__(
            prefix=["ClinVarResult-Set", "VariationArchive", "InterpretedRecord", "ClinicalAssertionList", "ClinicalAssertion"]
        )

    @parser_path(on_start=True)
    def new_record(self, elem):
        self.reset()
        self.latest.submitter_date = ClinVarXmlParser.parse_xml_date(elem.get("SubmissionDate"))

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
        "ReviewStatus")
    def parse_review_status(self, elem):
        review_status = elem.text
        self.latest.review_status = review_status
        self.latest.stars = CLINVAR_REVIEW_STATUS_TO_STARS.get(review_status, 0)

    @parser_path("Interpretation")
    def parse_date_last_evaluated(self, elem):
        if dle := elem.get("DateLastEvaluated"):
            self.latest.date_last_evaluated = ClinVarXmlParser.parse_xml_date(dle)

    @parser_path(
        "Interpretation",
        "Comment")
    def parse_interpretation_summary(self, elem):
        self.latest.interpretation_summary = elem.text

    @parser_path(
        "Interpretation",
        "Description")
    def parse_clinical_significance(self, elem):
        cs = elem.text
        self.latest.clinical_significance = CLINVAR_TO_VG_CLIN_SIG.get(cs, cs)

    @parser_path(
        "SimpleAllele",
        "OtherNameList",
        "Name"
    )
    def parse_possible_c_hgvs(self, elem):
        if elem.get("Type") is None and (hgvs := elem.text):
            self.latest.c_hgvs = ClinVarXmlParser.parse_hgvs(hgvs)

    @parser_path(
        "SimpleAllele",
        "AttributeSet",
        PP("Attribute", Type="HGVS"))
    def parse_c_hgvs(self, elem):
        if not self.latest.c_hgvs:
            if hgvs := elem.text:
                self.latest.c_hgvs = ClinVarXmlParser.parse_hgvs(hgvs)

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
        self.latest.variant_coordinate = f"{chr}:{start} {ref}>{alt} ({assembly})"


    @parser_path(
        PP("TraitSet", Type="DrugResponse"),
        PP("Trait", Type="DrugResponse"),
        "Name",
        PP("ElementValue", Type="Preferred"))
    def parse_drug_response(self, elem):
        if not self.latest.condition:
            self.latest.condition = elem.text

    @parser_path(
        PP("TraitSet", Type="Disease"),
        PP("Trait", Type="Disease"),
        "XRef")
    def parse_condition(self, elem):
        if not self.latest.condition or ":" not in self.latest.condition:
            if final_value := ClinVarXmlParser.clean_term_id_from_elem(elem):
                self.latest.condition = final_value

    @parser_path(
        PP("TraitSet", Type="DrugResponse"),
        PP("Trait", Type="DrugResponse"),
        "Name",
        PP("ElementValue", Type="Preferred"))
    def parse_drug_response(self, elem):
        if not self.latest.condition:
            self.latest.condition = elem.text


    @parser_path(
        PP("TraitSet", Type="Disease"),
        PP("Trait", Type="Disease"),
        "Name",
        PP("ElementValue", Type="Preferred"))
    def parse_condition_name(self, elem):
        if not self.latest.condition:
            condition_text = elem.text
            if condition_text.lower() not in ("not specified", "not provided"):
                self.latest.condition = elem.text
