from dataclasses import dataclass, field
from typing import List, TypedDict, Mapping, TypeVar, Generic, Any

import bs4
from django.contrib.auth.models import User
from django.db import models
from django.db.models import PROTECT, CASCADE
from django.db.models.signals import post_save
from django.dispatch import receiver
from guardian.shortcuts import assign_perm
from lazy import lazy
from model_utils.models import TimeStampedModel

from classification.enums import ShareLevel, SpecialEKeys
from classification.models import ClassificationModification, EvidenceKeyMap, classification_post_publish_signal, \
    Classification, flag_types, EvidenceKey, MultiCondition
from classification.models.evidence_mixin import VCDbRefDict
from annotation.regexes import DbRegexes
from flags.models import flag_comment_action, Flag, FlagComment, FlagResolution, FlagStatus
from genes.hgvs import CHGVS
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from ontology.models import OntologyTerm, OntologyService
from snpdb.models import Allele, Lab, GenomeBuild


@dataclass(frozen=True)
class ClinVarMessage:

    severity: str
    text: str

    @property
    def bs(self) -> str:
        if self.severity == "error":
            return "danger"
        elif self.severity == "warning":
            return "warning"
        elif self.severity == "info":
            return "info"
        return "info"


@dataclass(frozen=True)
class ClinVarMessages:

    messages: List[str] = field(default_factory=list)

    @staticmethod
    def error(message: str):
        return ClinVarMessages([ClinVarMessage(severity="error", text=message)])

    @staticmethod
    def warning(message: str):
        return ClinVarMessages([ClinVarMessage(severity="warning", text=message)])

    @staticmethod
    def info(message: str):
        return ClinVarMessages([ClinVarMessage(severity="info", text=message)])

    def __add__(self, other) -> 'ClinVarMessages':
        if not other:
            return self
        return ClinVarMessages(self.messages + other.messages)

    def __bool__(self):
        return bool(self.messages)

    def __iter__(self):
        return iter(self.messages)


empty_messages = ClinVarMessages()


class ClinVarSubmissionNotes:

    def __init__(self):
        self.errors: List[str] = list()

    def add_error(self, text):
        self.errors.append(text)


JsonType = TypeVar('JsonType')


class ValidatedJson(Generic[JsonType]):

    def __init__(self, json_data: JsonType, messages: ClinVarMessages = empty_messages):
        self.json_data = json_data
        self.messages = messages

    def to_json(self) -> JsonType:
        return self.json_data

    @staticmethod
    def _traverse_messages(json_data) -> ClinVarMessages:
        messages = empty_messages
        if isinstance(json_data, list):
            for val in json_data:
                messages += ValidatedJson._traverse_messages(val)
        elif isinstance(json_data, dict):
            for val in json_data.values():
                messages += ValidatedJson._traverse_messages(val)
        elif isinstance(json_data, ValidatedJson):
            messages += json_data.messages
            messages += ValidatedJson._traverse_messages(json_data.json_data)
        return messages

    @lazy
    def all_messages(self) -> ClinVarMessages:
        return ValidatedJson._traverse_messages(self)

    def __bool__(self):
        return bool(self.json_data) or bool(self.messages)


class ClinVarEvidenceKey:

    def __init__(self, evidence_key: EvidenceKey, value_obj: Any):
        self.evidence_key = evidence_key
        self.values = []
        self.messages = empty_messages

        value = None
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
                                self.messages += ClinVarMessages.error(f"ADMIN: \"{self.evidence_key.pretty_label}\" value of \"{sub_value}\" doesn't have a ClinVar value")
                            found = True
                            break
                    if not found:
                        self.values.append(sub_value)
                        self.messages += ClinVarMessages.error(f"\"{self.evidence_key.pretty_label}\" value of \"{sub_value}\" isn't a valid value, can't translate to ClinVar")
            else:
                raise ValueError(f"ADMIN: Trying to extract value from \"{self.evidence_key.pretty_label}\" that isn't a SELECT or MULTISELECT")

    def value(self, single: bool = True, optional: bool = False) -> ValidatedJson:
        messages = self.messages
        if len(self.values) == 0:
            if not optional:
                messages += ClinVarMessages.error(f"No value for required field \"{self.evidence_key.pretty_label}\"")
            return ValidatedJson(None if single else [], messages)
        elif len(self.values) == 1:
            return ValidatedJson(self.values[0] if single else self.values, messages)
        else:
            messages = self.messages
            if single:
                messages += ClinVarMessages.error(f"\"{self.evidence_key.pretty_label}\" expected single value, got multiple values")
            return ValidatedJson(self.values, messages)

    def __bool__(self):
        return len(self.values) > 0


class ClinVarCitation(TypedDict, total=False):
    db: str
    id: str
    url: str


class ClinVarExportStatus(models.TextChoices):
    SUBMIT_WHEN_READY = 'S', 'Submit When Ready'
    PENDING = 'P', 'Pending'
    REJECT = 'R', 'Reject'


class ClinVarExport(TimeStampedModel, GuardianPermissionsMixin):
    allele = models.ForeignKey(Allele, on_delete=PROTECT)
    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    classification_based_on = models.ForeignKey(ClassificationModification, on_delete=PROTECT)
    transcript = models.TextField()

    review_date = models.DateTimeField(null=True, blank=True)
    review_status = models.CharField(max_length=1, choices=ClinVarExportStatus.choices, default=ClinVarExportStatus.PENDING)
    dirty_date = models.DateTimeField(null=True, blank=True)

    submit_when_possible = models.BooleanField(default=False)
    withdrawn = models.BooleanField(default=False)

    CITATION_DB_MAPPING = {
        DbRegexes.PUBMED.db: "PubMed",
        DbRegexes.PMC.db: "pmc",
        DbRegexes.NCBIBookShelf.db: "BookShelf"  #TODO confirm this is correct mapping
    }

    @staticmethod
    def best_clinvar_candidate(cm1: ClassificationModification, cm2: ClassificationModification) -> ClassificationModification:
        # FIXME created date isn't the best, last curated would be more accurate
        if cm1.classification.created > cm2.classification.created:
            return cm1
        return cm2

    def update_with(self, cm: ClassificationModification) -> bool:
        is_update = self.id is None or self.classification_based_on != cm
        self.classification_based_on = cm
        # self.dirty_date = datetime.now()
        # only update dirty date if json that we will send changes
        return is_update

    @staticmethod
    def chgvs_for(cm: ClassificationModification):
        # FIXME, make 38 only
        if chgvs_str := cm.classification.chgvs_grch37:
            return CHGVS(chgvs_str)
        if chgvs_str := cm.get(SpecialEKeys.C_HGVS):
            return CHGVS(chgvs_str)
        return None

    @staticmethod
    def sync_allele(allele: Allele):

        new_count = 0
        updated_count = 0
        orphan_count = 0

        transcript_lab_to_cm = dict()

        cm: ClassificationModification
        # TODO - do we reduce this to just records shared globally??!
        # or do we accept both, and treat the share level as a to be confirmed globally
        for cm in ClassificationModification.objects.filter(
                classification__withdrawn=False,
                classification__variant__variantallele__allele=allele,
                is_last_published=True,
                share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).select_related('classification'):

            classification = cm.classification
            # decision, grabbing the 38 representation
            if c_parts := ClinVarExport.chgvs_for(cm):
                transcript_no_version = c_parts.transcript_parts.identifier
                lab_str = str(classification.lab_id)
                db_key = f"{lab_str}*{transcript_no_version}"
                use_cm = cm
                if existing := transcript_lab_to_cm.get(db_key):
                    use_cm = ClinVarExport.best_clinvar_candidate(cm, existing)
                transcript_lab_to_cm[db_key] = use_cm

        cve: ClinVarExport
        existing_records = {f"{cve.lab.id}*{cve.transcript}": cve for cve in ClinVarExport.objects.filter(allele=allele).select_related('lab')}
        for lab_transcript, cm in transcript_lab_to_cm.items():
            cve: ClinVarExport
            if cve := existing_records.pop(lab_transcript, None):
                if cve.update_with(cm) or cve.withdrawn:
                    cve.withdrawn = False
                    cve.save()
                    updated_count += 1
            else:
                c_parts = ClinVarExport.chgvs_for(cm)
                cve = ClinVarExport(
                    allele=allele,
                    lab=cm.classification.lab,
                    transcript=lab_transcript.split("*", maxsplit=1)[1]
                )
                cve.update_with(cm)
                cve.save()
                new_count += 1

        for orphan in existing_records.values():
            if not orphan.withdrawn:
                # TODO allow user to withdraw
                orphan.withdrawn = True
                orphan.save()
                orphan_count += 1

        return {
            "new": new_count,
            "updated": updated_count,
            "withdrawn": orphan_count
        }

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
        pubmed_refs = {ref.get('id'): ref for ref in self.classification_based_on.db_refs if ref.get('db') in ClinVarExport.CITATION_DB_MAPPING}
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
                soup = bs4.BeautifulSoup(value, 'lxml')
                return soup.text
        return value

    @lazy
    def as_json(self) -> ValidatedJson:
        data = dict()
        data["assertionCriteria"] = self.json_assertion_criteria
        data["clinicalSignificance"] = self.json_clinical_significance
        data["conditionSet"] = self.condition_set
        data["localId"] = ValidatedJson("TODO", ClinVarMessages.error("ADMIN: localID & lockKey are under development"))
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
                return ValidatedJson(None, ClinVarMessages.error(f"No normalised c.hgvs in genome build {genome_build}"))
        except:
            return ValidatedJson(None, ClinVarMessages.error(f"Could not determine genome build of submission"))

    @property
    def json_record_status(self) -> ValidatedJson:
        return ValidatedJson("novel", ClinVarMessages.warning("ADMIN: Defaulting status to novel, need to store if previously submitted"))

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
            data["citation"] = ValidatedJson(acmg_citation, ClinVarMessages.warning(f"AssertionMethod value \"{assertion_criteria}\", providing default ACMG as citation"))
        if method := EvidenceKeyMap.cached_key(SpecialEKeys.ASSERTION_METHOD).pretty_value(assertion_criteria):
            data["method"] = method
        else:
            data["method"] = ValidatedJson(None, ClinVarMessages.error("No value provided for Assertion Method"))

        return ValidatedJson(data)

    @property
    def json_clinical_significance(self) -> ValidatedJson:
        data = dict()
        if citations := self.citation_refs:
            def citation_to_json(citation: VCDbRefDict) -> ClinVarCitation:
                return {
                    "db": ClinVarExport.CITATION_DB_MAPPING.get(citation.get("db")),
                    "id": str(citation.get("idx"))  # TODO confirm this is the kind of ID they want, not the prefixed one?
                }
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
        messages = ClinVarMessages()
        condition_list = list()
        data["condition"] = condition_list
        if conditions := self.classification_based_on.classification.condition_resolution_obj:
            if len(conditions.terms) >= 2:
                if join := conditions.join:
                    join = MultiCondition(join)
                    if join != MultiCondition.CO_OCCURRING:
                        messages += ClinVarMessages.error("ClinVar API only supports Co-occurring for multiple messages")

            def condition_to_json(condition: OntologyTerm) -> ValidatedJson:
                # supported "OMIM", "MedGen", "Orphanet", "MeSH", "HP", "MONDO"
                messages = empty_messages
                if condition.ontology_service not in (OntologyService.OMIM, OntologyService.ORPHANET, OntologyService.HPO, OntologyService.MONDO):
                    messages += ClinVarMessages.error(f"Ontology \"{condition.ontology_service}\" is not supported by ClinVar")
                return ValidatedJson({
                    "db": condition.ontology_service,
                    "id": condition.id
                }, messages)

            for condition in conditions.terms:
                condition_list.append(condition_to_json(condition))
        else:
            messages += ClinVarMessages.error("No standard condition terms")
        return ValidatedJson(data, messages)

    @property
    def observed_in(self) -> ValidatedJson:
        data = dict()
        data["affectedStatus"] = self.clinvar_value(SpecialEKeys.AFFECTED_STATUS).value(single=True)
        if allele_origin := self.clinvar_value(SpecialEKeys.ALLELE_ORIGIN).value(single=True, optional=True):
            data["allele_origin"] = allele_origin
        else:
            data["alleleOrigin"] = ValidatedJson("germline", ClinVarMessages.info("Defaulting \"AlleleOrigin\" to germline as no value provided"))
        data["collectionMethod"] = "curation"  # TODO confirm hardcoded of curation
        # numberOfIndividuals do we do anything with this?
        return ValidatedJson(data)


@receiver(post_save, sender=ClinVarExport)
def set_condition_alias_permissions(sender, created: bool, instance: ClinVarExport, **kwargs):  # pylint: disable=unused-argument
    if created:
        group = instance.lab.group
        assign_perm(ClinVarExport.get_read_perm(), group, instance)
        assign_perm(ClinVarExport.get_write_perm(), group, instance)


class ClinVarExportSubmission(TimeStampedModel, GuardianPermissionsMixin):
    clinvar_export = models.ForeignKey(ClinVarExport, on_delete=CASCADE)
    evidence = models.JSONField()
    submission_status = models.TextField()


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender,
              classification: Classification,
              previously_published: ClassificationModification,
              newly_published: ClassificationModification,
              previous_share_level: ShareLevel,
              user: User,
              **kwargs):

    cve: ClinVarExport
    if cve := ClinVarExport.objects.filter(classification_based_on__classification=classification).first():
        cve.update_with(newly_published)
        cve.save()


@receiver(flag_comment_action, sender=Flag)
def check_for_discordance(sender, flag_comment: FlagComment, old_resolution: FlagResolution, **kwargs):  # pylint: disable=unused-argument
    """
    Keeps condition_text_match in sync with the classifications when withdraws happen/finish
    """
    flag = flag_comment.flag
    if flag.flag_type == flag_types.classification_flag_types.classification_withdrawn:
        cl: Classification
        if cl := Classification.objects.filter(flag_collection=flag.collection.id).first():
            cve: ClinVarExport
            if cve := ClinVarExport.objects.filter(classification_based_on__classification=cl).first():
                cve.withdrawn = flag_comment.resolution.status == FlagStatus.OPEN
                cve.save()
