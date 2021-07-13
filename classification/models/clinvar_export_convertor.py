from typing import List, Any, Mapping, TypedDict
import bs4
from lazy import lazy

from annotation.regexes import DbRegexes
from classification.enums import SpecialEKeys
from classification.json_utils import JSON_MESSAGES_EMPTY, JsonMessages, ValidatedJson
from classification.models import ClassificationModification, EvidenceKeyMap, EvidenceKey, \
    MultiCondition
from classification.models.evidence_mixin import VCDbRefDict
from ontology.models import OntologyTerm, OntologyService
from snpdb.models import GenomeBuild


class ClinVarEvidenceKey:
    """
    Handles the basic mapping for evidence keys to ClinVar
    Revolves around select/mutli-select fields having equivilents
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


class ClinVarCitation(TypedDict, total=False):
    db: str
    id: str
    url: str


class ClinVarExportConverter:

    CITATION_DB_MAPPING = {
        DbRegexes.PUBMED.db: "PubMed",
        DbRegexes.PMC.db: "pmc",
        DbRegexes.NCBIBookShelf.db: "BookShelf"  # TODO confirm this is correct mapping
    }

    def __init__(self, classification_based_on: ClassificationModification):
        self.classification_based_on = classification_based_on

    # convenience method for forms
    def evidence_key(self, key: str):
        return EvidenceKeyMap.cached().get(key).pretty_value(self.classification_based_on.get(key), dash_for_none=True)

    @property
    def allele_origin(self):
        return self.evidence_key(SpecialEKeys.ALLELE_ORIGIN)

    @property
    def mode_of_inheritance(self):
        return self.evidence_key(SpecialEKeys.MODE_OF_INHERITANCE)

    @property
    def affected_status(self):
        return self.evidence_key(SpecialEKeys.AFFECTED_STATUS)

    @property
    def genome_build(self):
        return self.evidence_key(SpecialEKeys.GENOME_BUILD)

    @property
    def interpretation_summary(self):
        # strip out any XML
        interpret = self.evidence_key(SpecialEKeys.INTERPRETATION_SUMMARY)
        soup = bs4.BeautifulSoup(interpret, 'lxml')
        return soup.text

    @property
    def curation_context(self):
        return self.evidence_key(SpecialEKeys.CURATION_CONTEXT)

    @property
    def assertion_method(self):
        return self.evidence_key(SpecialEKeys.ASSERTION_METHOD)

    # probably wont use this
    @property
    def patient_phenotype(self):
        return self.evidence_key(SpecialEKeys.PATIENT_PHENOTYPE)

    @property
    def clinical_significance(self):
        return self.evidence_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)

    @property
    def citation_refs(self) -> List[VCDbRefDict]:
        # TODO allow other database references "PubMed", "BookShelf", "DOI", "pmc"
        # TODO check to see if the PubMeds exist?
        pubmed_refs = {ref.get('id'): ref for ref in self.classification_based_on.db_refs if ref.get('db') in ClinVarExportConverter.CITATION_DB_MAPPING}
        unique_refs = list(pubmed_refs.values())
        unique_refs.sort(key=lambda x: x.get('id'))
        return unique_refs

    @property
    def condition(self) -> str:
        """
        Note condition isn't actually sent, but the value from ConditionTextMatching
        this is just here for reference
        """
        return self.evidence_key(SpecialEKeys.CONDITION)

    @property
    def curated_date(self):
        return self.classification_based_on.curated_date

    @property
    def c_hgvs(self):
        gb = GenomeBuild.get_from_fuzzy_string(self.genome_build)
        return self.classification_based_on.classification.get_c_hgvs(genome_build=gb)

    def clinvar_value(self, key: str) -> ClinVarEvidenceKey:
        evidence_key = EvidenceKeyMap.cached_key(key)
        value = self.classification_based_on.evidence.get(key)
        return ClinVarEvidenceKey(evidence_key, value)

    def value(self, key: str) -> Any:
        if value := self.classification_based_on.get(key):
            if isinstance(value, str):
                # FIXME, only do if key is multiline text
                soup = bs4.BeautifulSoup(value, 'lxml')
                return soup.text
        return value

    @lazy
    def as_json(self) -> ValidatedJson:
        data = dict()
        data["assertionCriteria"] = self.json_assertion_criteria
        data["clinicalSignificance"] = self.json_clinical_significance
        data["conditionSet"] = self.condition_set
        data["localId"] = ValidatedJson("TODO", JsonMessages.error("ADMIN: localID & lockKey are under development"))
        data["observedIn"] = self.observed_in
        data["recordStatus"] = self.json_record_status
        data["releaseStatus"] = "public"
        data["variantSet"] = self.variant_set

        return ValidatedJson(data)

    @property
    def variant_set(self) -> ValidatedJson:
        try:
            genome_build = self.classification_based_on.get_genome_build()
            if c_hgvs := self.classification_based_on.classification.get_c_hgvs(genome_build):
                return ValidatedJson({"hgvs": c_hgvs})
            else:
                return ValidatedJson(None, JsonMessages.error(f"No normalised c.hgvs in genome build {genome_build}"))
        except:
            return ValidatedJson(None, JsonMessages.error("Could not determine genome build of submission"))

    @property
    def json_record_status(self) -> ValidatedJson:
        return ValidatedJson("novel", JsonMessages.warning("ADMIN: Defaulting status to novel, need to store if previously submitted"))

    @property
    def json_assertion_criteria(self) -> ValidatedJson:
        data = dict()
        acmg_citation = {
            "db": "PubMed",
            "id": "PubMed:25741868"
        }
        assertion_criteria = self.value(SpecialEKeys.ASSERTION_METHOD)
        if assertion_criteria == "acmg":
            data["citation"] = acmg_citation
        else:
            data["citation"] = ValidatedJson(acmg_citation, JsonMessages.warning(f"AssertionMethod value \"{assertion_criteria}\", providing default ACMG as citation"))
        if method := EvidenceKeyMap.cached_key(SpecialEKeys.ASSERTION_METHOD).pretty_value(assertion_criteria):
            data["method"] = method
        else:
            data["method"] = ValidatedJson(None, JsonMessages.error("No value provided for Assertion Method"))

        return ValidatedJson(data)

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
            data["citations"] = [citation_to_json(citation) for citation in citations]
        data["clinicalSignificanceDescription"] = self.clinvar_value(SpecialEKeys.CLINICAL_SIGNIFICANCE).value(single=True)
        if interpret := self.value(SpecialEKeys.INTERPRETATION_SUMMARY):
            data["comment"] = interpret
        if date_last_evaluated := self.value(SpecialEKeys.CURATION_DATE):
            # FIXME also check validation date?
            data["dateLastEvaluated"] = date_last_evaluated
        if mode_of_inheritance := self.clinvar_value(SpecialEKeys.MODE_OF_INHERITANCE):
            data["modeOfInheritance"] = mode_of_inheritance.value(single=False)
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
        data = dict()
        data["affectedStatus"] = self.clinvar_value(SpecialEKeys.AFFECTED_STATUS).value(single=True)
        if allele_origin := self.clinvar_value(SpecialEKeys.ALLELE_ORIGIN).value(single=True, optional=True):
            data["allele_origin"] = allele_origin
        else:
            data["alleleOrigin"] = ValidatedJson("germline", JsonMessages.info("Defaulting \"AlleleOrigin\" to germline as no value provided"))
        data["collectionMethod"] = "curation"  # TODO confirm hardcoded of curation
        # numberOfIndividuals do we do anything with this?
        return ValidatedJson(data)