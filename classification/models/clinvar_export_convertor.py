from re import RegexFlag
from typing import List, Any, Mapping, TypedDict, Union
import bs4
from lazy import lazy

from annotation.regexes import DbRegexes
from classification.enums import SpecialEKeys, EvidenceKeyValueType
from genes.hgvs import CHGVS
from library.utils import html_to_text
from uicore.json.validated_json import JsonMessages, JSON_MESSAGES_EMPTY, ValidatedJson
from classification.models import ClassificationModification, EvidenceKeyMap, EvidenceKey, \
    MultiCondition, ClinVarExport, classification_flag_types
from classification.models.evidence_mixin import VCDbRefDict
from ontology.models import OntologyTerm, OntologyService
from snpdb.models import ClinVarKey
import re
import json


# Code in this file is responsible for converting VariantGrid formatted classifications to ClinVar JSON


class ClinVarEvidenceKey:
    """
    Handles the basic mapping of a single field to ClinVar
    Revolves around select/multi-select fields having equivalents stored in "clinvar" attribute of the EvidenceKey
    """

    def __init__(self, evidence_key: EvidenceKey, value_obj: Any):
        self.evidence_key = evidence_key
        self.values = []
        self.messages = JSON_MESSAGES_EMPTY

        value: Any
        if isinstance(value_obj, Mapping):
            value = value_obj.get('value')
        else:
            value = value_obj
        if not value:
            value = []
        elif not isinstance(value, list):
            value = [value]

        if values := value:
            if options := evidence_key.virtual_options:
                for sub_value in values:
                    found = False
                    for option in options:
                        if option.get('key') == sub_value:
                            if clinvar := option.get('clinvar'):
                                self.values.append(clinvar)
                            else:
                                self.values.append(sub_value)
                                self.messages += JsonMessages.error(f"ADMIN: \"{self.evidence_key.pretty_label}\" value of \"{sub_value}\" doesn't have a ClinVar value")
                            found = True
                            break
                    if not found:
                        self.values.append(sub_value)
                        self.messages += JsonMessages.error(f"\"{self.evidence_key.pretty_label}\" value of \"{sub_value}\" isn't a valid value, can't translate to ClinVar")
            else:
                raise ValueError(f"ADMIN: Trying to extract value from \"{self.evidence_key.pretty_label}\" that isn't a SELECT or MULTISELECT")

    def value(self, single: bool = True, optional: bool = False) -> ValidatedJson:
        """
        :param single: Is a single value expected. Adds an error to ValidatedJson if number of values is 2 or more.
        :param optional: Is the value optional, if not adds an error to ValidatedJson if number of values is zero.
        """
        messages = self.messages
        if len(self.values) == 0:
            if not optional:
                messages += JsonMessages.error(f"No value for required field \"{self.evidence_key.pretty_label}\"")
            return ValidatedJson(None if single else [], messages)
        elif len(self.values) == 1:
            return ValidatedJson(self.values[0] if single else self.values, messages)
        else:
            messages = self.messages
            if single:
                messages += JsonMessages.error(f"\"{self.evidence_key.pretty_label}\" expected single value, got multiple values")
            return ValidatedJson(self.values, messages)

    def __bool__(self):
        return len(self.values) > 0


# Dictionary definitions, we don't have many since we deal more with ValidatedJSon where a typed dictionary doesn't fit


class ClinVarCitation(TypedDict, total=False):
    """
    Dictionary representation of a citation in the ClinVar format
    """
    db: str
    id: str
    url: str


class ClinVarExportConverter:

    CITATION_DB_MAPPING = {
        DbRegexes.PUBMED.db: "PubMed",
        DbRegexes.PMC.db: "pmc",
        DbRegexes.NCBIBookShelf.db: "BookShelf"
    }

    FLAG_TYPES_TO_MESSAGES = {
        classification_flag_types.classification_withdrawn: JsonMessages.error("Classification has since been withdrawn"),
        classification_flag_types.transcript_version_change_flag: JsonMessages.error("Classification has open transcript version flag"),
        classification_flag_types.matching_variant_warning_flag: JsonMessages.error("Classification has open variant warning flag"),
        classification_flag_types.discordant: JsonMessages.error("Classification is in discordance"),
        classification_flag_types.internal_review: JsonMessages.error("Classification is in internal review"),
    }

    def __init__(self, clinvar_export_record: ClinVarExport):
        """
        :param clinvar_export_record: the record to use as a basis to convert to ClinVar json (with validation)
        """
        self.clinvar_export_record = clinvar_export_record

    @lazy
    def clinvar_key(self) -> ClinVarKey:
        return self.clinvar_export_record.clinvar_allele.clinvar_key

    @property
    def classification_based_on(self) -> ClassificationModification:
        return self.clinvar_export_record.classification_based_on

    @property
    def citation_refs(self) -> List[VCDbRefDict]:
        pubmed_refs = {ref.get('id'): ref for ref in self.classification_based_on.db_refs if ref.get('db') in ClinVarExportConverter.CITATION_DB_MAPPING}
        unique_refs = list(pubmed_refs.values())
        unique_refs.sort(key=lambda x: x.get('id'))
        return unique_refs

    def clinvar_value(self, key: str) -> ClinVarEvidenceKey:
        evidence_key = EvidenceKeyMap.cached_key(key)
        value = self.classification_based_on.evidence.get(key)
        return ClinVarEvidenceKey(evidence_key, value)

    def value(self, key: str) -> Any:
        value = self.classification_based_on.get(key)
        if isinstance(value, str) and EvidenceKeyMap.cached_key(key).value_type == EvidenceKeyValueType.TEXT_AREA:
            return html_to_text(value, preserve_lines=True)
        else:
            return value

    @lazy
    def as_validated_json(self) -> ValidatedJson:
        data = dict()
        if self.classification_based_on is None:
            return ValidatedJson(None, JsonMessages.error("No classification is currently associated with this allele and condition"))
        else:
            data["assertionCriteria"] = self.json_assertion_criteria
            data["clinicalSignificance"] = self.json_clinical_significance
            data["conditionSet"] = self.condition_set
            allele_id = self.clinvar_export_record.clinvar_allele.allele_id

            local_id = f"ALLELE_{allele_id}"
            c = self.classification_based_on.classification
            local_key = c.lab.group_name + "/" + c.lab_record_id

            """
            # Below was used when condition was literally part of the localID
            # Important, using the umbrella term as it doesn't change
            condition = self.clinvar_export_record.condition_resolved
            term_ids = "_".join([term.id.replace(":", "_") for term in condition.terms])
            if join := condition.join:
                join = MultiCondition(join)
                term_ids = f"{term_ids}_{join.value}"
            """

            data["localID"] = local_id
            data["localKey"] = local_key
            data["observedIn"] = self.observed_in
            data["releaseStatus"] = "public"
            data["variantSet"] = self.variant_set

            messages = JSON_MESSAGES_EMPTY
            for flag in self.classification_based_on.classification.flag_collection.flags(only_open=True):
                if message := ClinVarExportConverter.FLAG_TYPES_TO_MESSAGES.get(flag.flag_type):
                    messages += message

            return ValidatedJson(data, messages)

    @property
    def variant_set(self) -> ValidatedJson:
        try:
            genome_build = self.classification_based_on.get_genome_build()
            if c_hgvs := self.classification_based_on.classification.get_c_hgvs(genome_build):
                c_hgvs_obj = CHGVS(c_hgvs)
                c_hgvs_no_gene = c_hgvs_obj.without_gene_symbol_str

                variant_data = [{"hgvs": c_hgvs_no_gene}]
                json_data = {"variant": [{"hgvs": c_hgvs_no_gene}]}

                hgvs_errors = JSON_MESSAGES_EMPTY
                if c_hgvs.startswith("ENST"):
                    hgvs_errors += JsonMessages.error("ClinVar doesn't support Ensembl transcripts")

                gene_symbols = list()
                if gene_symbol := self.value(SpecialEKeys.GENE_SYMBOL):
                    gene_symbols.append({"symbol": gene_symbol})

                json_data = {
                    "variant": [
                        {
                            "hgvs": ValidatedJson(c_hgvs_no_gene,hgvs_errors),
                            "gene": gene_symbols
                        }
                    ]
                }

                return json_data
            else:
                return ValidatedJson(None, JsonMessages.error(f"No normalised c.hgvs in genome build {genome_build}"))
        except BaseException:
            return ValidatedJson(None, JsonMessages.error("Could not determine genome build of submission"))

    @property
    def json_assertion_criteria(self) -> Union[dict, ValidatedJson]:
        assertion_criteria = self.value(SpecialEKeys.ASSERTION_METHOD)

        assertion_method_lookups = self.clinvar_key.assertion_method_lookup
        acmg_criteria = {
            "citation": {
                "db": "PubMed",
                "id": "PubMed:25741868"
            },
            "method": EvidenceKeyMap.cached_key(SpecialEKeys.ASSERTION_METHOD).pretty_value("acmg")
        }

        if assertion_criteria == "acmg":
            return acmg_criteria
        else:
            for key, criteria in assertion_method_lookups.items():
                raw_criteria = criteria
                if criteria == "acmg":
                    criteria = acmg_criteria

                if not assertion_criteria:
                    if not criteria or re.compile(key, RegexFlag.IGNORECASE).match(""):
                        return ValidatedJson(criteria, JsonMessages.info(f"Using config for assertion method \"{key}\" : {json.dumps(raw_criteria)}"))
                else:
                    expr = re.compile(key, RegexFlag.IGNORECASE)
                    if expr.match(assertion_criteria):
                        return ValidatedJson(criteria, JsonMessages.info(f"Using config for assertion method \"{key}\" : {json.dumps(raw_criteria)}"))
            else:
                return ValidatedJson(None, JsonMessages.error(f"No match for assertion method of \"{assertion_criteria}\""))

    @property
    def json_clinical_significance(self) -> ValidatedJson:
        data = dict()
        if citations := self.citation_refs:
            def citation_to_json(citation: VCDbRefDict) -> ClinVarCitation:
                citation: ClinVarCitation = {
                    "db": ClinVarExportConverter.CITATION_DB_MAPPING.get(citation.get("db")),
                    "id": str(citation.get("idx"))  # TODO confirm this is the kind of ID they want, not the prefixed one?
                }
                return citation
            data["citation"] = [citation_to_json(citation) for citation in citations]
        data["clinicalSignificanceDescription"] = self.clinvar_value(SpecialEKeys.CLINICAL_SIGNIFICANCE).value(single=True)

        comment_parts: List[str] = list()

        if interpret := self.value(SpecialEKeys.INTERPRETATION_SUMMARY):
            comment_parts.append(interpret.strip())

        if self.clinvar_key.inject_acmg_description and (acmg_summary := self.classification_based_on.criteria_strength_summary()):
            comment_parts.append(acmg_summary)

        if comment_parts:
            data["comment"] = "\n\n".join(comment_parts)

        if date_last_evaluated := self.value(SpecialEKeys.CURATION_DATE):
            # FIXME fall back to other date types? or at least raising a warning
            data["dateLastEvaluated"] = date_last_evaluated

        if mode_of_inheritance := self.clinvar_value(SpecialEKeys.MODE_OF_INHERITANCE):
            data["modeOfInheritance"] = mode_of_inheritance.value(single=True)
        return ValidatedJson(data)

    @property
    def condition_set(self) -> ValidatedJson:
        data = dict()
        messages = JSON_MESSAGES_EMPTY
        condition_list = list()
        data["condition"] = condition_list
        if conditions := self.classification_based_on.classification.condition_resolution_obj:
            if len(conditions.terms) >= 2:
                if join := conditions.join:
                    join = MultiCondition(join)
                    if join != MultiCondition.CO_OCCURRING:
                        messages += JsonMessages.error("ClinVar API only supports Co-occurring for multiple messages")

            def condition_to_json(condition: OntologyTerm) -> ValidatedJson:
                # supported "OMIM", "MedGen", "Orphanet", "MeSH", "HP", "MONDO"
                messages = JSON_MESSAGES_EMPTY
                if condition.ontology_service not in (OntologyService.OMIM, OntologyService.ORPHANET, OntologyService.HPO, OntologyService.MONDO):
                    messages += JsonMessages.error(f"Ontology \"{condition.ontology_service}\" is not supported by ClinVar")
                return ValidatedJson({
                    "db": condition.ontology_service,
                    "id": condition.id
                }, messages)

            for condition in conditions.terms:
                condition_list.append(condition_to_json(condition))
        else:
            messages += JsonMessages.error("No standard condition terms")
        return ValidatedJson(data, messages)

    @property
    def observed_in(self) -> ValidatedJson:
        # can return array, but we only have one
        # (though
        data = dict()
        affected_status_value = self.clinvar_value(SpecialEKeys.AFFECTED_STATUS)
        if affected_status_value:
            data["affectedStatus"] = affected_status_value.value(single=True)
        elif default_affected_status := self.clinvar_key.default_affected_status:
            data["affectedStatus"] = ValidatedJson(default_affected_status, JsonMessages.info("Using configured default for affected status"))
        else:
            data["affectedStatus"] = ValidatedJson(None, JsonMessages.error("Affected status not provided and no default config provided"))
        if allele_origin := self.clinvar_value(SpecialEKeys.ALLELE_ORIGIN).value(single=True, optional=True):
            data["alleleOrigin"] = allele_origin
        else:
            data["alleleOrigin"] = ValidatedJson("germline", JsonMessages.info("Defaulting \"Allele origin\" to \"germline\" as no value provided"))
        data["collectionMethod"] = "clinical testing"
        # numberOfIndividuals do we do anything with this?
        return ValidatedJson([data])
