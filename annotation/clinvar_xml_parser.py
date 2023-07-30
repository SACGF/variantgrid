from dataclasses import dataclass
from typing import Optional
import re

from library.utils import ExportRow, export_column
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

@dataclass
class ConvertedClinVarOutput(ExportRow):
    record_id: Optional[str] = None
    genome_build: Optional[str] = None
    c_hgvs: Optional[str] = None
    condition: Optional[str] = None
    clinical_significance: Optional[str] = None
    gene_symbol: Optional[str] = None
    date_last_evaluated: Optional[str] = None
    interpretation_summary: Optional[str] = None
    assertion_method: Optional[str] = None
    allele_origin: Optional[str] = None

    @export_column("lab_record_id")
    def export_record_id(self):
        return self.record_id

    @export_column("genome_build")
    def genome_build(self):
        return self.genome_build

    @export_column("c_hgvs")
    def export_c_hgvs(self):
        return self.c_hgvs

    @export_column("gene_symbol")
    def export_gene_symbol(self):
        return self.gene_symbol

    @export_column("condition")
    def export_condition(self):
        return self.condition or "Not set"

    @export_column("clinical_significance")
    def export_clinical_significance(self):
        return CLINVAR_TO_VG_CLIN_SIG.get(self.clinical_significance, "VUS")

    # @export_column("allele_origin")
    # def export_allele_origin(self):
    #     return self.allele_origin or "unknown"

    @export_column("interpretation_summary")
    def export_interpretation_summary(self):
        return self.interpretation_summary

    @property
    def is_good_quality(self) -> bool:
        # if not (self.c_hgvs and self.gene_symbol):
        #     return False
        # if "?" in self.c_hgvs or "[" in self.c_hgvs or " " in self.c_hgvs:
        #     return False
        return True

    def __repr__(self):
        return str(self.to_csv())


class ClinVarParser(XmlParser):

    RE_GOOD_CHGVS = re.compile("^(N._[0-9]+[.][0-9]+:c[.][0-9_a-zA-Z>]+)( .*)?$")
    RE_ORPHA = re.compile("ORPHA([0-9]+)")

    def __init__(self, genome_build: str):
        self.latest: Optional[ConvertedClinVarOutput] = None
        self.genome_build = genome_build
        super().__init__()

    def reset(self):
        if self.latest and self.latest.is_good_quality:
            self.set_yieldable(self.latest)
        self.latest = ConvertedClinVarOutput()
        self.latest.genome_build = self.genome_build

    def finish(self):
        self.reset()

    @parser_path("ClinVarResult-Set", "ClinVarSet", "ClinVarAssertion")
    def new_record(self, elem):
        self.reset()

    @parser_path("ClinVarResult-Set", "ClinVarSet", "ClinVarAssertion", PP("ClinVarAccession", Type="SCV"))
    def record_id(self, elem):
        self.latest.record_id = elem.get("Acc")

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
            if match := ClinVarParser.RE_GOOD_CHGVS.match(hgvs):
                hgvs = match.group(1)
                self.latest.c_hgvs = hgvs

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
                if m := ClinVarParser.RE_ORPHA.match(id):
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
