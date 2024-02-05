#
# As ClinVar XML layout has changed since this was first written, commenting it out as it's out of date
#
from typing import List

# import json
#
# from Bio import Entrez
#
# from annotation.clinvar_xml_parser import ClinVarXmlParser, ClinVarXmlParserOutput, CLINVAR_REVIEW_STATUS_TO_STARS, \
#     CLINVAR_TO_VG_CLIN_SIG
# from annotation.models import ClinVarRecord
# from library.utils.xml_utils import parser_path, PP
#
#
# class ClinVarXmlParserViaRCVs(ClinVarXmlParser):
#     """
#     This implementation is deprecated, as it requires us to make multiple requests to get the same data as via VCV
#     Left this code in case this is helpeful for reading clinvar xml downloads etc
#     """
#
#     PARSER_VERSION = 100
#
#     @classmethod
#     def load_from_clinvar_id(cls, clinvar_variation_id: int) -> ClinVarXmlParserOutput:
#         cv_handle = Entrez.esummary(db="clinvar", retmode="json", id=clinvar_variation_id)
#         json_data = json.loads(cv_handle.read())
#         cv_handle.close()
#
#         all_urls: list[str] = [f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={clinvar_variation_id}&retmode=json"]
#         all_rcvs: list[str] = []
#         parsed_results: list[ClinVarRecord] = []
#         if result := json_data.get("result"):
#             if uuids := result.get("uids"):
#                 if uuid_data := result.get(uuids[0]):
#                     if supporting_submissions := uuid_data.get("supporting_submissions"):
#                         if rcvs := supporting_submissions.get("rcv"):
#                             # rcv is a combination of a condition and variant of which there's 1 or more
#                             # we need to retrieve all of them to get the full data
#                             all_rcvs = rcvs
#         if all_rcvs:
#             all_urls += [f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&id={rcv}" for rcv in all_rcvs]
#             handle = Entrez.efetch(db="clinvar", rettype="clinvarset", id=all_rcvs)
#             parsed_results = cls.load_from_input(handle)
#             handle.close()
#
#         return ClinVarXmlParserOutput(
#             urls=all_urls,
#             all_records=parsed_results
#         )
#
#     def __init__(self):
#         super().__init__(
#             prefix=["ClinVarResult-Set", "ClinVarSet", "ClinVarAssertion"]
#         )
#
#     @parser_path(on_start=True)
#     def new_record(self, elem):
#         self.reset()
#
#     @parser_path(PP("ClinVarAccession", Type="SCV"))
#     def record_id(self, elem):
#         self.latest.record_id = elem.get("Acc")
#         self.latest.org_id = elem.get("OrgID")
#         self.latest.date_clinvar_created = ClinVarXmlParser.parse_xml_date(elem.get("DateCreated"))
#         self.latest.date_clinvar_updated = ClinVarXmlParser.parse_xml_date(elem.get("DateUpdated"))
#
#     @parser_path(
#         "ObservedIn",
#         "Sample",
#         "Origin")
#     def parse_allele_origin(self, elem):
#         self.latest.allele_origin = elem.text
#
#     @parser_path(
#         PP("MeasureSet", Type="Variant"),
#         PP("Measure", Type="Variation"),
#         "AttributeSet",
#         PP("Attribute", Type="HGVS"))
#     def parse_c_hgvs(self, elem):
#         self.assign_better_hgvs(elem.text)
#
#     @parser_path(
#         PP("MeasureSet", Type="Variant"),
#         PP("Measure", Type="Variation"),
#         "SequenceLocation")
#     def parse_variant_coordinate(self, elem):
#         assembly = elem.get("Assembly")
#         chro = elem.get("Chr")
#         start = elem.get("start")
#         ref = elem.get("referenceAllele")
#         alt = elem.get("alternateAllele")
#         self.latest.variant_coordinate = f"{chro}:{start} {ref}>{alt} ({assembly})"
#
#     @parser_path(
#         "ClinVarSubmissionID")
#     def parse_reviewer(self, elem):
#         self.latest.submitter = elem.get("submitter")
#         if submitter_date := elem.get("submitterDate"):
#             self.latest.submitter_date = ClinVarXmlParser.parse_xml_date(submitter_date)
#
#     @parser_path(
#         "ClinicalSignificance",
#         "ReviewStatus")
#     def parse_review_status(self, elem):
#         review_status = elem.text
#         self.latest.review_status = review_status
#         self.latest.stars = CLINVAR_REVIEW_STATUS_TO_STARS.get(review_status, 0)
#
#     @parser_path(
#         PP("MeasureSet", Type="Variant"),
#         PP("Measure", Type="Variation"),
#         PP("MeasureRelationship", Type="variant in gene"),
#         "Symbol",
#         PP("ElementValue", Type="Preferred"))
#     def parse_gene_symbol(self, elem):
#         self.latest.gene_symbol = elem.text
#
#     @parser_path(
#         "ClinicalSignificance",
#         "Description")
#     def parse_clinical_significance(self, elem):
#         cs = elem.text
#         if cs := cs.lower():
#             self.latest.clinical_significance = CLINVAR_TO_VG_CLIN_SIG.get(cs, cs)
#
#     @parser_path(
#         "ClinicalSignificance",
#         "Comment")
#     def parse_interpretation_summary(self, elem):
#         self.latest.interpretation_summary = elem.text
#
#     @parser_path(
#         "AttributeSet",
#         PP("Attribute", Type="AssertionMethod"))
#     def parse_assertion_method(self, elem):
#         self.latest.assertion_method = elem.text
#
#     @parser_path(
#         "ClinicalSignificance")
#     def parse_date_last_evaluated(self, elem):
#         if last_evaluated := elem.get("DateLastEvaluated"):
#             self.latest.date_last_evaluated = ClinVarXmlParser.parse_xml_date(last_evaluated)
#
#     @parser_path(
#         PP("TraitSet", Type="DrugResponse"),
#         PP("Trait", Type="DrugResponse"),
#         "Name",
#         PP("ElementValue", Type="Preferred"))
#     def parse_drug_response(self, elem):
#         if not self.latest.condition:
#             self.latest.condition = elem.text
#
#     @parser_path(
#         PP("TraitSet", Type="Disease"),
#         PP("Trait", Type="Disease"),
#         "XRef")
#     def parse_condition(self, elem):
#         if not self.latest.condition:
#             if final_value := ClinVarXmlParser.clean_term_id_from_elem(elem):
#                 self.latest.condition = final_value
