import json
import re
from re import RegexFlag
from typing import List, Any, Mapping, TypedDict, Union

from lazy import lazy

from annotation.regexes import DbRegexes
from classification.enums import SpecialEKeys, EvidenceKeyValueType, ShareLevel
from classification.models import ClassificationModification, EvidenceKeyMap, EvidenceKey, \
    MultiCondition, ClinVarExport, classification_flag_types, Classification
from classification.models.evidence_mixin import VCDbRefDict
from genes.hgvs import CHGVS
from library.utils import html_to_text
from ontology.models import OntologyTerm, OntologyService
from snpdb.models import ClinVarKey
from uicore.json.validated_json import JsonMessages, JSON_MESSAGES_EMPTY, ValidatedJson

# Code in this file is responsible for converting VariantGrid formatted classifications to ClinVar JSON
CLINVAR_ACCEPTED_TRANSCRIPTS = {"NM_", "NR_"}

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

    def __len__(self):
        return len(self.values)


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
        classification_flag_types.classification_outstanding_edits: JsonMessages.error("Classification has un-submitted changes"),
        classification_flag_types.classification_pending_changes: JsonMessages.error("Classification has pending changes"),
        # classification_flag_types.classification_not_public - this has special handling to include a comment
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

    @staticmethod
    def citation_to_json(citation: VCDbRefDict) -> ClinVarCitation:
        db = ClinVarExportConverter.CITATION_DB_MAPPING.get(citation.get("db"))
        id_part = citation['id'].replace(' ', '')
        if db == "PubMed":
            # Special support for PubMed to be PMID
            # (at some point should fix that in the citation JSON)
            id_part = f"PMID:{citation['idx']}"

        citation: ClinVarCitation = {
            "db": db,
            "id": id_part
        }
        return citation

    @staticmethod
    def condition_to_json(condition: OntologyTerm) -> ValidatedJson:
        # supported "OMIM", "MedGen", "Orphanet", "MeSH", "HP", "MONDO"
        messages = JSON_MESSAGES_EMPTY
        if condition.ontology_service not in (OntologyService.OMIM, OntologyService.ORPHANET,
                                              OntologyService.HPO, OntologyService.MONDO):
            messages += JsonMessages.error(f"Ontology \"{condition.ontology_service}\" is not supported by ClinVar")

        """
        # Examples
        OMIM and 100800
        MeSH and D000130
        Orphanet and ORPHA155
        MedGen and C0001080
        Mondo and MONDO:0015263
        """

        id_part = condition.id
        if condition.ontology_service == OntologyService.OMIM:
            id_part = str(condition.index)  # OMIM is not 0 prefixed
        elif condition.ontology_service == OntologyService.ORPHANET:
            id_part = f"ORPHA{str(condition.index)}"  # ORPHA is not 0 prefixed

        return ValidatedJson({
            "db": condition.ontology_service,
            "id": id_part
        }, messages)

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

            data["localID"] = local_id
            data["localKey"] = local_key
            data["observedIn"] = self.observed_in
            data["releaseStatus"] = "public"
            data["variantSet"] = self.variant_set

            messages = JSON_MESSAGES_EMPTY
            for flag in self.classification_based_on.classification.flag_collection.flags(only_open=True):
                if flag.flag_type == classification_flag_types.classification_not_public:
                    message_text = "Not submitting to ClinVar as record is flagged to be excluded"
                    if last_comment := flag.flagcomment_set.order_by('-created').first():
                        if last_comment_text := last_comment.text:
                            if "Not submitting to ClinVar as the record matched the following pattern(s):" in last_comment_text:
                                message_text = last_comment_text
                            else:
                                message_text += " - " + last_comment_text
                    messages += JsonMessages.error(message_text)

                elif message := ClinVarExportConverter.FLAG_TYPES_TO_MESSAGES.get(flag.flag_type):
                    messages += message

            # see if other shared classifications for the clinvar_key variant combo don't have a resolved condition
            # but only if they don't have an open don't share flag
            allele = self.clinvar_export_record.clinvar_allele.allele
            if other_classifications_for_key := Classification.objects.filter(
                    withdrawn=False,
                    variant__in=allele.variants,
                    share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS,
                    lab__clinvar_key=self.clinvar_key
            ).exclude(id=self.classification_based_on.id):
                for c in other_classifications_for_key:
                    has_condition = (resolved_condition := c.condition_resolution_obj) and len(resolved_condition.terms) >= 1
                    if not has_condition:
                        # make sure it doesn't have an exclude flag, no point complaining about that
                        if not c.flag_collection_safe.get_open_flag_of_type(flag_type=classification_flag_types.classification_not_public):
                            messages += JsonMessages.error(f"Another classification for this allele '{c.lab_record_id}' has an unresolved condition with text '{c.get(SpecialEKeys.CONDITION)}'")

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
                for accepted_transcript in CLINVAR_ACCEPTED_TRANSCRIPTS:
                    if c_hgvs.startswith(accepted_transcript):
                        break
                else:
                    hgvs_errors += JsonMessages.error(f"ClinVar only accepts transcripts starting with one of {CLINVAR_ACCEPTED_TRANSCRIPTS}")

                gene_symbols = list()
                if gene_symbol := self.value(SpecialEKeys.GENE_SYMBOL):
                    gene_symbols.append({"symbol": gene_symbol})

                json_data = {
                    "variant": [
                        {
                            "hgvs": ValidatedJson(c_hgvs_no_gene, hgvs_errors),
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
                "id": "PMID:25741868"
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
            data["citation"] = [ClinVarExportConverter.citation_to_json(citation) for citation in citations]
        data["clinicalSignificanceDescription"] = self.clinvar_value(SpecialEKeys.CLINICAL_SIGNIFICANCE).value(single=True)

        comment_parts: List[str] = list()

        if interpret := self.value(SpecialEKeys.INTERPRETATION_SUMMARY):
            comment_parts.append(interpret.strip())

        if self.clinvar_key.inject_acmg_description and (acmg_summary := self.classification_based_on.criteria_strength_summary()):
            comment_parts.append(f"ACMG/AMP criteria applied: {acmg_summary}")

        if comment_parts:
            full_comment = " ".join(comment_parts)
            full_comment = full_comment.replace("\n", " ")
            data["comment"] = full_comment

        if date_last_evaluated := self.value(SpecialEKeys.CURATION_DATE):
            data["dateLastEvaluated"] = date_last_evaluated
        elif date_last_reviewed := self.value(SpecialEKeys.CURATION_VERIFIED_DATE):
            data["dateLastEvaluated"] = ValidatedJson(date_last_reviewed, JsonMessages.warning("No curation date provided, falling back to curation verified date."))

        messages = JsonMessages()
        if mode_of_inheritance := self.clinvar_value(SpecialEKeys.MODE_OF_INHERITANCE):
            # TODO might need the concept of NO_KEY so we can have "mode_of_inheritance": NO_KEY
            # so we can still put the warning against mode_of_inheritance, but not actually produce it in the real JSON
            # as ClinVar doesn't accept modeOfInheritance: null
            if len(mode_of_inheritance) > 1:
                messages += JsonMessages.warning("ClinVar only accepts a single value for mode of inheritance. There are multiple values for mode of inheritance against this record, so omitting this field.")
            else:
                data["modeOfInheritance"] = mode_of_inheritance.value(single=True)
        return ValidatedJson(data, messages)

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

            for condition in conditions.terms:
                condition_list.append(ClinVarExportConverter.condition_to_json(condition))
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
