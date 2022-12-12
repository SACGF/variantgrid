from dataclasses import dataclass
from typing import List, Any, Mapping, TypedDict, Optional
from lazy import lazy
from annotation.models import CitationFetchRequest
from annotation.models.models_citations import CitationSource
from annotation.regexes import db_citation_regexes
from classification.enums import SpecialEKeys, EvidenceKeyValueType, ShareLevel
from classification.models import ClassificationModification, EvidenceKeyMap, EvidenceKey, \
    MultiCondition, ClinVarExport, classification_flag_types, Classification, ClinVarExportStatus, \
    ClinVarExportSubmission, CLINVAR_EXPORT_CONVERSION_VERSION
from genes.hgvs import CHGVS
from library.utils import html_to_text, JsonObjType, JsonDiffs
from ontology.models import OntologyTerm, OntologyService, OntologyTermStatus
from snpdb.models import ClinVarKey, ClinVarCitationsModes
from uicore.json.validated_json import JsonMessages, JSON_MESSAGES_EMPTY, ValidatedJson

# Code in this file is responsible for converting VariantGrid formatted classifications to ClinVar JSON
CLINVAR_ACCEPTED_TRANSCRIPTS = {"NM_", "NR_"}


class ClinVarCitation(TypedDict, total=False):
    """
    Dictionary representation of a citation in the ClinVar format
    """
    db: str
    id: str
    url: str


@dataclass(frozen=True)
class ClinVarExportChanges:
    """
    Keeps track of what changed between what was last submitted and current.
    Consider renaming to ClinVarExportDifferences
    """

    grouping_changes: Optional[JsonDiffs]
    """ Changes in the grouping header, None if the previous record didn't have a grouping """

    body_changes: Optional[JsonDiffs]
    """ Changes in the body """

    version_change: bool = False
    """ Indicates if the record previously submitted under a different version """

    @property
    def grouping_changes_json(self):
        return self.grouping_changes.to_json("previous", "current")

    @property
    def body_changes_json(self):
        return self.body_changes.to_json("previous", "current")

    def __bool__(self):
        return self.version_change or bool(self.grouping_changes) or bool(self.body_changes)


@dataclass(frozen=True)
class ClinVarExportData:
    """
    Provides functionality around a ClinVarExport to update it or put it into a batch
    """
    clinvar_export: 'ClinVarExport'
    """ The record this export data is for """

    grouping: Optional[ValidatedJson] = None
    """ The body data we'll be using (could either be as it's cached in the record or what the most up to date grouping
    header would look like"""

    body: Optional[ValidatedJson] = None
    """ The body data we'll be using (could either be as it's cached in the record or what the most up to date body
    would look like"""

    @property
    def no_record(self) -> bool:
        """ Indicates if there's no valid classification """
        return self.clinvar_export.classification_based_on is None

    @lazy
    def body_json(self) -> JsonObjType:
        return self.body.pure_json()

    @property
    def has_errors(self) -> bool:
        if self.grouping and self.grouping.has_errors:
            return True
        if self.body and self.body.has_errors:
            return True
        return False

    @lazy
    def _previous_submission(self) -> ClinVarExportSubmission:
        return self.clinvar_export.previous_submission

    @property
    def is_new(self):
        """ Indicates if this is the first time the record has been submitted.
            Important, this is not the same as already having an SCV """
        return not (self.clinvar_export.pk and self._previous_submission)

    @property
    def is_valid(self):
        """ Are we able to send this data to ClinVar - important does NOT check status """
        return bool(self.grouping) and bool(self.body) and not self.has_errors and 'assertionCriteria' in self.grouping

    @lazy
    def changes(self) -> Optional[ClinVarExportChanges]:
        """
        Differences between this data and the previously submitted
        """
        if self.no_record:
            return None
        if previous_submission := self._previous_submission:
            grouping_changes = JsonDiffs.differences(previous_submission.submission_grouping or dict(), self.grouping.pure_json())
            body_changes = JsonDiffs.differences(previous_submission.submission_body, self.body.pure_json())

            return ClinVarExportChanges(
                grouping_changes=grouping_changes,
                body_changes=body_changes,
                version_change=previous_submission.submission_version != CLINVAR_EXPORT_CONVERSION_VERSION
            )
        return None

    def submission_full(self) -> JsonObjType:
        """
        JSON data including recordStatus and scv, not normally included in body as we don't want the assigning of an SCV
        to generate a change
        :return:
        """
        json_data = self.body.pure_json()
        if scv := self.clinvar_export.scv:
            json_data["recordStatus"] = "update"
            json_data["clinvarAccession"] = scv
        else:
            json_data["recordStatus"] = "novel"
        return json_data

    @property
    def assertion_criteria(self) -> Optional[ClinVarCitation]:
        return self.grouping.pure_json().get("assertionCriteria") if self.grouping else None

    @property
    def local_id(self) -> str:
        return self.body_json.get("localID")

    @property
    def local_key(self) -> str:
        return self.body_json.get("localKey")

    def apply(self):
        """
        Updates the underlying ClinVarExport record with this data, and updates the status
        """
        clinvar_export = self.clinvar_export

        status: ClinVarExportStatus
        if (cm := clinvar_export.classification_based_on) and cm.classification.flag_collection_safe.get_open_flag_of_type(
            classification_flag_types.classification_not_public):
            status = ClinVarExportStatus.EXCLUDE
        elif self.has_errors:
            status = ClinVarExportStatus.IN_ERROR
        else:
            if self.is_new:
                status = ClinVarExportStatus.NEW_SUBMISSION
            elif self.changes:
                status = ClinVarExportStatus.CHANGES_PENDING
            else:
                status = ClinVarExportStatus.UP_TO_DATE

        clinvar_export.status = status
        clinvar_export.submission_grouping_validated = self.grouping.serialize()
        clinvar_export.submission_body_validated = self.body.serialize()
        lazy.invalidate(clinvar_export, 'submission_grouping')
        lazy.invalidate(clinvar_export, 'submission_body')
        clinvar_export.save()

    @staticmethod
    def current(clinvar_export: ClinVarExport) -> 'ClinVarExportData':
        """
        Make ClinVarExoprtData based on current evidence (good if we just reloaded clinvar_export from the DB)
        """
        return ClinVarExportData(
            clinvar_export=clinvar_export,
            grouping=clinvar_export.submission_grouping,
            body=clinvar_export.submission_body
        )


class ClinVarEvidenceKey:
    """
    Handles the basic mapping of a single field to ClinVar
    Revolves around select/multi-select fields having equivalents stored in "clinvar" attribute of the EvidenceKey
    """

    def __init__(self, evidence_key: EvidenceKey, value_obj: Any):
        self.evidence_key = evidence_key
        self.valid_values = []
        self.invalid_values = []
        self.conversion_messages = JSON_MESSAGES_EMPTY

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
                                self.valid_values.append(clinvar)
                            else:
                                self.invalid_values.append(sub_value)
                                self.conversion_messages += JsonMessages.warning(f"\"{self.evidence_key.pretty_label}\" value of \"{sub_value}\" doesn't have a ClinVar equivalent.")
                            found = True
                            break
                    if not found:
                        self.invalid_values.append(sub_value)
                        self.conversion_messages += JsonMessages.warning(f"\"{self.evidence_key.pretty_label}\" value of \"{sub_value}\" isn't valid, can't translate to ClinVar.")
            else:
                raise ValueError(f"ADMIN: Trying to extract value from \"{self.evidence_key.pretty_label}\" that isn't a SELECT or MULTISELECT.")

    def value(self, single: bool = True, optional: bool = False) -> ValidatedJson:
        """
        :param single: Is a single value expected. Adds an error to ValidatedJson if number of values is 2 or more.
        :param optional: Is the value optional, if not adds an error to ValidatedJson if number of values is zero.
        """
        messages = self.conversion_messages
        if single and len(self) > 1:
            messages += JsonMessages.warning(f"\"{self.evidence_key.pretty_label}\" within ClinVar only accepts a single value. Record has the values {self.valid_values + self.invalid_values} for this field.")

        if not optional:
            if messages:
                # if mandatory, then any warnings about conversion have to be upgraded to errors
                messages += JsonMessages.error(f"\"{self.evidence_key.pretty_label}\" is mandatory but has conversion issues.")
            elif not self.valid_values:
                # if mandatory, we obviously need a value
                messages += JsonMessages.error(f"\"{self.evidence_key.pretty_label}\" is mandatory but has no value.")

        if messages:
            # mandatory or not, if we've had conversion warnings, provide no value to be safe (even if there's potentially valid values)
            return ValidatedJson.make_void(messages)
        elif len(self) == 0:
            # even when a key is optional, ClinVar sometimes doesn't like the value null being provided, and instead prefers the key to be omitted
            return ValidatedJson.make_void()
        else:
            return ValidatedJson(self.valid_values[0] if single else self.valid_values)

    def __len__(self):
        return len(self.valid_values) + len(self.invalid_values)

    def __bool__(self):
        # think this is default behaviour for bool
        return len(self) > 0


class ClinVarExportConverter:

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

    @staticmethod
    def is_exclude_citation(citation_json: ClinVarCitation) -> bool:
        if pmid_id := citation_json.get('id'):
            return pmid_id in {
                "PMID:25741868",  # ACMG criteria
                "PMID:28492532",  # Sherloc criteria
                "PMID:30192042",  # Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion
            }

    @property
    def clinvar_key(self) -> ClinVarKey:
        return self.clinvar_export_record.clinvar_allele.clinvar_key

    @property
    def classification_based_on(self) -> ClassificationModification:
        return self.clinvar_export_record.classification_based_on

    CITATION_SOURCE_TO_CLINVAR = {
        CitationSource.PUBMED: "PubMed",
        CitationSource.PUBMED_CENTRAL: "pmc",
        CitationSource.NCBI_BOOKSHELF: "Bookshelf"  # FIXME, don't know if this is what NCBI wants, waiting on email
    }

    @property
    def citations(self) -> List[ValidatedJson]:
        request_ids = list()
        if self.clinvar_key.citations_mode == ClinVarCitationsModes.interpretation_summary_only:
            if text := self.classification_based_on.get(SpecialEKeys.INTERPRETATION_SUMMARY):
                request_ids = db_citation_regexes.search(text)
        else:
            # non citation refs will be ignored
            request_ids = self.classification_based_on.db_refs

        citation_list = list()
        for citation in CitationFetchRequest.fetch_all_now(request_ids).all_citations:
            if clinvar_db := ClinVarExportConverter.CITATION_SOURCE_TO_CLINVAR.get(citation.source):
                # NEED to test for NBCI
                citation_json = {
                    "db": clinvar_db,
                    "id": citation.id
                }
                messages = JSON_MESSAGES_EMPTY

                if ClinVarExportConverter.is_exclude_citation(citation_json):
                    citation_list.append(ValidatedJson.make_void(JsonMessages.info(f"{citation_json.get('id')} is omitted as it's not specific for this variant.")))
                else:
                    if citation.error:
                        messages += JsonMessages.error(f"Could not retrieve \"{citation.id}\", might not be valid.")

                    citation_list.append(ValidatedJson(citation_json, messages=messages))

        return citation_list

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
    def condition_to_json(condition: OntologyTerm) -> ValidatedJson:
        # supported "OMIM", "MedGen", "Orphanet", "MeSH", "HP", "MONDO"
        messages = JSON_MESSAGES_EMPTY
        if condition.ontology_service not in (OntologyService.OMIM, OntologyService.ORPHANET,
                                              OntologyService.HPO, OntologyService.MONDO):
            messages += JsonMessages.error(f"Ontology \"{condition.ontology_service}\" is not supported by ClinVar")

        if condition.ontology_service in (OntologyService.OMIM, OntologyService.HPO, OntologyService.MONDO) and not condition.is_valid_for_condition:
            # we don't keep Orphanet records, so just hope for the best with them
            
            if condition.status == OntologyTermStatus.STUB:
                messages += JsonMessages.error(f"We have no record of condition \"{condition}\". If this is valid it will be automatically resolved when our system refreshes ontology terms next.")
            elif condition.status == OntologyTermStatus.DEPRECATED:
                messages += JsonMessages.error(
                    f"Term \"{condition}\" has been marked as deprecated. You will need to condition match to a newer term to submit this data to ClinVar.")
            elif condition.status == OntologyTermStatus.NON_CONDITION:
                messages += JsonMessages.error(
                    f"Term \"{condition}\" is not valid to use as a condition under curation. You will need to condition match to a newer term to submit this data to ClinVar."
                )

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

    def convert(self) -> ClinVarExportData:

        if self.classification_based_on is None:
            return ClinVarExportData(
                clinvar_export=self.clinvar_export_record
            )
        else:
            grouping = ValidatedJson({"assertionCriteria": self.json_assertion_criteria})

            data = {}
            data["clinicalSignificance"] = self.json_clinical_significance
            data["conditionSet"] = self.condition_set
            allele_id = self.clinvar_export_record.clinvar_allele.allele_id

            local_id = f"ALLELE_{allele_id}"
            c = self.classification_based_on.classification
            local_key = c.lab.group_name + "/" + c.lab_record_id

            data["localID"] = local_id
            data["localKey"] = local_key
            data["observedIn"] = self.observed_in
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

            return ClinVarExportData(
                clinvar_export=self.clinvar_export_record,
                grouping=grouping,
                body=ValidatedJson(data, messages)
            )

    @property
    def variant_set(self) -> ValidatedJson:
        try:
            genome_build = self.classification_based_on.get_genome_build()
            if c_hgvs := self.classification_based_on.classification.get_c_hgvs(genome_build):
                c_hgvs_obj = CHGVS(c_hgvs)
                c_hgvs_no_gene = c_hgvs_obj.without_gene_symbol_str

                hgvs_errors = JSON_MESSAGES_EMPTY
                for accepted_transcript in CLINVAR_ACCEPTED_TRANSCRIPTS:
                    if c_hgvs.startswith(accepted_transcript):
                        break
                else:
                    hgvs_errors += JsonMessages.error(f"ClinVar only accepts transcripts starting with one of {CLINVAR_ACCEPTED_TRANSCRIPTS}")

                gene_symbols = []
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

                return ValidatedJson(json_data)
            else:
                return ValidatedJson(None, JsonMessages.error(f"No normalised c.hgvs in genome build {genome_build}"))
        except BaseException:
            return ValidatedJson(None, JsonMessages.error("Could not determine genome build of submission"))

    @property
    def json_assertion_criteria(self) -> ValidatedJson:
        assertion_criteria = self.value(SpecialEKeys.ASSERTION_METHOD) or ""
        if mapped_assertion_method := self.clinvar_key.assertion_criteria_vg_to_code(assertion_criteria):
            if "url" in mapped_assertion_method or ("db" in mapped_assertion_method and "id" in mapped_assertion_method):
                return ValidatedJson(mapped_assertion_method, JsonMessages.info(f"Mapped from \"{assertion_criteria}\"."))
            else:
                return ValidatedJson(mapped_assertion_method, JsonMessages.error("ADMIN: Mapped to invalid assertionCriteria"))
        else:
            return ValidatedJson(mapped_assertion_method or "", JsonMessages.error(f'Could not map \"{assertion_criteria}\" to an assertionCriteria.'))

    @property
    def json_clinical_significance(self) -> ValidatedJson:
        data = {}
        if citations := self.citations:
            data["citation"] = citations
        data["clinicalSignificanceDescription"] = self.clinvar_value(SpecialEKeys.CLINICAL_SIGNIFICANCE).value(single=True)

        comment_parts: List[str] = []

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
            data["modeOfInheritance"] = mode_of_inheritance.value(single=True, optional=True)

        return ValidatedJson(data, messages)

    @property
    def condition_set(self) -> ValidatedJson:
        data = {}
        messages = JSON_MESSAGES_EMPTY
        condition_list = []
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
        data = {}
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

    @staticmethod
    def clinvar_export_data(clinvar_export: ClinVarExport, update: bool) -> ClinVarExportData:
        if update:
            data = ClinVarExportConverter(clinvar_export).convert()
            data.apply()
        # always load fresh just to be safe
        return ClinVarExportData.current(clinvar_export)