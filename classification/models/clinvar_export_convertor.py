import re
from dataclasses import dataclass
from functools import cached_property
from typing import Any, Mapping, TypedDict, Optional

from annotation.models import CitationFetchRequest
from annotation.models.models_citations import CitationSource
from annotation.regexes import db_citation_regexes
from classification.enums import SpecialEKeys, EvidenceKeyValueType, ShareLevel, AlleleOriginBucket
from classification.models import ClassificationModification, EvidenceKeyMap, EvidenceKey, \
    MultiCondition, ClinVarExport, classification_flag_types, Classification, ClinVarExportStatus, \
    ClinVarExportSubmission, CLINVAR_EXPORT_CONVERSION_VERSION, CuratedDate
from genes.hgvs import CHGVS
from library.utils import html_to_text, JsonObjType, JsonDiffs, invalidate_cached_property
from ontology.models import OntologyTerm, OntologyService, OntologyTermStatus
from snpdb.models import ClinVarKey, ClinVarCitationsModes, ClinVarExportTypeBucket, ClinVarExportAssertionMethod, \
    ClinVarCitationSource
from uicore.json.validated_json import JsonMessages, JSON_MESSAGES_EMPTY, ValidatedJson

"""
Code in this file is responsible for converting VariantGrid formatted classifications to ClinVar JSON
It also applies data level validation
"""

# ClinVar doesn't accept ENST for example
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

    grouping: ValidatedJson
    """ JSON that has to match between records to be put into the same clinvar batch """

    body: ValidatedJson
    """ The body data we'll be using (could either be as it's cached in the record or what the most up to date body
    would look like"""

    @property
    def no_record(self) -> bool:
        """ Indicates if there's no valid classification """
        return self.clinvar_export.classification_based_on is None

    @cached_property
    def body_json(self) -> JsonObjType:
        return self.body.pure_json()

    @property
    def has_errors(self) -> bool:
        if self.grouping and self.grouping.has_errors:
            return True
        if self.body and self.body.has_errors:
            return True
        return False

    @cached_property
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
        if self.clinvar_export.delete_pending and self.clinvar_export.scv:
            return True
        return self.assertion_criteria and bool(self.body) and not self.has_errors

    @cached_property
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
        if not self.clinvar_export.delete_pending:
            if scv := self.clinvar_export.scv:
                json_data["recordStatus"] = "update"
                json_data["clinvarAccession"] = scv
            else:
                json_data["recordStatus"] = "novel"
        return json_data

    @property
    def assertion_criteria(self) -> Optional[ClinVarCitation]:
        if grouping := self.grouping:
            if pure_json := grouping.pure_json():
                return pure_json.get("assertionCriteria")
        return None

    @property
    def local_id(self) -> str:
        return self.body_json.get("localID")

    @property
    def local_key(self) -> str:
        return self.body_json.get("local_key")

    def apply(self):
        """
        Updates the underlying ClinVarExport record with this data and updates the status
        """
        clinvar_export = self.clinvar_export

        status: ClinVarExportStatus
        if self.clinvar_export.deleted:
            status = ClinVarExportStatus.DELETED
        elif self.clinvar_export.delete_pending:
            if self.changes:
                status = ClinVarExportStatus.CHANGES_PENDING
            else:
                status = ClinVarExportStatus.UP_TO_DATE
        elif (cm := clinvar_export.classification_based_on) and cm.classification.flag_collection_safe.get_open_flag_of_type(
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
        invalidate_cached_property(clinvar_export, 'submission_grouping')
        invalidate_cached_property(clinvar_export, 'submission_body')
        clinvar_export.save()

    @staticmethod
    def current(clinvar_export: ClinVarExport) -> 'ClinVarExportData':
        """
        Make ClinVarExportData based on current evidence (good if we just reloaded clinvar_export from the DB)
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
        # these flags are now handled by AlleleInfo
        # classification_flag_types.transcript_version_change_flag: JsonMessages.error("Classification has open transcript version flag"),
        # classification_flag_types.matching_variant_warning_flag: JsonMessages.error("Classification has open variant warning flag"),
        classification_flag_types.discordant: JsonMessages.error("Classification is in discordance"),
        classification_flag_types.internal_review: JsonMessages.error("Classification is in internal review"),
        classification_flag_types.classification_outstanding_edits: JsonMessages.error("Classification has un-submitted changes"),
        classification_flag_types.classification_pending_changes: JsonMessages.error("Classification has pending changes"),
        # classification_flag_types.classification_not_public - this has special handling to include a comment
    }

    BOOKSHELF_ID_RE = re.compile(".*(NBK[0-9]+)")

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
        return False

    @property
    def clinvar_key(self) -> ClinVarKey:
        return self.clinvar_export_record.clinvar_allele.clinvar_key

    @property
    def classification_based_on(self) -> ClassificationModification:
        return self.clinvar_export_record.classification_based_on

    CITATION_SOURCE_TO_CLINVAR = {
        CitationSource.PUBMED: ClinVarCitationSource.PUBMED,
        CitationSource.PUBMED_CENTRAL: ClinVarCitationSource.PUBMED_CENTRAL,
        CitationSource.NCBI_BOOKSHELF: ClinVarCitationSource.NCBI_BOOKSHELF
    }

    @property
    def citations(self) -> list[ValidatedJson]:
        request_ids = []
        if self.clinvar_key.citations_mode == ClinVarCitationsModes.interpretation_summary_only:
            if text := self.classification_based_on.get(SpecialEKeys.INTERPRETATION_SUMMARY):
                request_ids = db_citation_regexes.search(text)
        else:
            # non-citation refs will be ignored
            request_ids = self.classification_based_on.db_refs

        citation_list = []
        for citation in CitationFetchRequest.fetch_all_now(request_ids).all_citations:
            if clinvar_db := ClinVarExportConverter.CITATION_SOURCE_TO_CLINVAR.get(citation.source):
                citation_id = citation.id
                if citation.source == CitationSource.NCBI_BOOKSHELF:
                    if m := ClinVarExportConverter.BOOKSHELF_ID_RE.match(citation_id):
                        citation_id = m.group(1)

                citation_json = {
                    "db": clinvar_db.value,
                    "id": citation_id
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

    def convert_date(self, value: str) -> ValidatedJson:
        # doesn't so much convert a date, as makes sure it's YYYY-MM-DD and is valid, e.g., not 2024-13-43
        if not CuratedDate.convert_classification_date_str(value):
            return ValidatedJson(value, JsonMessages.error("Invalid date"))
        return ValidatedJson(value)

    @staticmethod
    def condition_to_json(condition: OntologyTerm) -> ValidatedJson:
        # supported "OMIM", "MedGen", "Orphanet", "MeSH", "HP", "MONDO"
        messages = JSON_MESSAGES_EMPTY
        if condition.ontology_service not in (OntologyService.OMIM, OntologyService.ORPHANET,
                                              OntologyService.HPO, OntologyService.MONDO):
            messages += JsonMessages.error(f"Ontology \"{condition.ontology_service}\" is not supported by ClinVar")

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
        db = condition.ontology_service
        if condition.ontology_service == OntologyService.OMIM:
            id_part = str(condition.index)  # OMIM is not 0 prefixed
        elif condition.ontology_service == OntologyService.ORPHANET:
            id_part = f"ORPHA{str(condition.index)}"  # ORPHA is not 0 prefixed
            db = "Orphanet"

        return ValidatedJson({
            "db": db,
            "id": id_part
        }, messages)

    def convert(self) -> ClinVarExportData:

        if self.clinvar_export_record.deleted:
            return ClinVarExportData(
                clinvar_export=self.clinvar_export_record,
                body=ValidatedJson.make_void(JsonMessages.info("This record has been deleted from ClinVar - no further syncing is possible")),
                grouping=ValidatedJson.make_void()
            )
        elif self.clinvar_export_record.delete_pending:
            data = {}
            # note we inject accession into non-deleted records live, and accession into deleted records straight into the DB
            # this is because this data is used to calculate differences from any previous submission, and we don't want

            if scv := self.clinvar_export_record.scv:
                data["accession"] = scv
            else:
                data["accession"] = ValidatedJson.make_void(JsonMessages.error("To delete, there has to be a SCV"))
            if reason := self.clinvar_export_record.delete_reason:
                data["reason"] = reason

            return ClinVarExportData(
                clinvar_export=self.clinvar_export_record,
                body=ValidatedJson(data),
                grouping=ValidatedJson.make_void()
            )
        if self.classification_based_on is None:
            return ClinVarExportData(
                clinvar_export=self.clinvar_export_record,
                body=ValidatedJson(None, messages=JsonMessages.error("No classification is currently associated with this allele and condition")),
                grouping=ValidatedJson({})
            )
        else:
            # remove warnings as if two records share the same assertionCriteria but have differently worded
            # warnings, that shouldn't affect anything
            grouping = ValidatedJson({
                "assertionCriteria": self.json_assertion_criteria,
            }).without_warnings()

            allele_id = self.clinvar_export_record.clinvar_allele.allele_id

            local_id = f"ALLELE_{allele_id}"
            c = self.classification_based_on.classification
            local_key = f"CR_{c.id}"

            # note clinVarAssertion gets injected later
            # as the SCV can change outside of the JSON generation
            # and we don't want to detect a change in JSON when the only change is SCV has being assigned
            data = {
                "conditionSet": self.condition_set,
                "localID": local_id,
                "local_key": local_key,
                "observedIn": self.observed_in,
                "variantSet": self.variant_set
            }
            match self.clinvar_export_record.clinvar_allele.clinvar_export_bucket:
                case ClinVarExportTypeBucket.GERMLINE:
                    data["germlineClassification"] = self.germline_classification
                case ClinVarExportTypeBucket.ONCOGENIC:
                    data["oncogenicityClassification"] = self.oncogenicity_classification
                case ClinVarExportTypeBucket.CLINICAL_IMPACT:
                    data["clinicalImpactClassification"] = self.clinical_impact_classification

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

            if not self.classification_based_on.classification.allele_info.latest_validation.include:
                messages += JsonMessages.error("There are outstanding variant matching warnings for this record")

            # see if other shared classifications for the clinvar_key variant combo don't have a resolved condition
            # but only if they don't have an open don't share flag
            allele = self.clinvar_export_record.clinvar_allele.allele

            if other_classifications_for_key := Classification.objects.filter(
                    withdrawn=False,
                    allele_info__allele=allele,
                    share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS,
                    allele_origin_bucket=self.classification_based_on.classification.allele_origin_bucket,
                    lab__clinvar_key=self.clinvar_key
            ).exclude(id=self.classification_based_on.id):
                for c in other_classifications_for_key:
                    has_condition = (resolved_condition := c.condition_resolution_obj) and len(resolved_condition.terms) >= 1
                    if not has_condition:
                        # make sure it doesn't have an exclude flag, no point complaining about that
                        if not c.flag_collection_safe.get_open_flag_of_type(flag_type=classification_flag_types.classification_not_public):
                            messages += JsonMessages.error(f'Another classification for this allele & allele origin "{c.cr_lab_id}" has an unresolved condition with text "{c.get(SpecialEKeys.CONDITION)}"')

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
                else:
                    gene_symbols = ValidatedJson([], JsonMessages.error("No gene symbol provided"))

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
    def mapped_assertion_method(self) -> Optional[tuple[str, ClinVarExportAssertionMethod]]:
        export_type_bucket = self.clinvar_export_record.clinvar_allele.clinvar_export_bucket
        assertion_criterias = list(sorted(self.classification_based_on.get_value_list(SpecialEKeys.ASSERTION_METHOD)))
        if len(assertion_criterias) == 0:
            assertion_criterias = [""]

        # pick the assertion method with the priority closest to 1
        assertion_method: Optional[str, ClinVarExportAssertionMethod] = None
        assertion_method_priority: Optional[int] = None
        for assertion_criteria in assertion_criterias:
            for mapping in self.clinvar_key.clinvarexportassertionmethodmapping_set.filter(
                    assertion_method__export_type=export_type_bucket).order_by('order'):
                if mapping.matches(assertion_criteria):
                    if assertion_method_priority is None or assertion_method_priority > mapping.order:
                        assertion_method = assertion_criteria, mapping.assertion_method
                        assertion_method_priority = mapping.order
                        break

        if assertion_method:
            return assertion_method
        else:
            for assertion_criteria in assertion_criterias:
                for possible_assertion_method in ClinVarExportAssertionMethod.objects.filter(export_type=export_type_bucket).order_by('id').all():
                    if possible_assertion_method.matches(assertion_criteria):
                        return assertion_criteria, possible_assertion_method
        return None

    @property
    def json_assertion_criteria(self) -> ValidatedJson:
        if mapped_assertion_method_tuple := self.mapped_assertion_method:
            return ValidatedJson(mapped_assertion_method_tuple[1].as_clinvar_json, messages=JsonMessages.info(f"Mapped assertion method \"{mapped_assertion_method_tuple[0]}\" to \"{mapped_assertion_method_tuple[1].label}\""))
        else:
            assertion_criterias = [x for x in sorted(self.classification_based_on.get_value_list(SpecialEKeys.ASSERTION_METHOD)) if x]
            assertion_criterias_str = ", ".join(f"\"{ac}\"" for ac in assertion_criterias) if assertion_criterias else "-blank-"
            export_type_bucket = ClinVarExportTypeBucket(self.clinvar_export_record.clinvar_allele.clinvar_export_bucket).name
            return ValidatedJson.make_void(JsonMessages.error(f"Could not map assertion methods {assertion_criterias_str} to {export_type_bucket} assertion method"))

    @property
    def base_classification(self) -> ValidatedJson:
        data = {}
        # we exclude citations unless we include interpretation summary
        if self.clinvar_key.include_interpretation_summary:
            if citations := self.citations:
                data["citation"] = citations

        comment_parts: list[str] = []
        if self.clinvar_key.include_interpretation_summary and (interpret := self.value(SpecialEKeys.INTERPRETATION_SUMMARY)):
            comment_parts.append(interpret.strip())

        if self.clinvar_key.inject_acmg_description and (criteria_summary := self.classification_based_on.criteria_strength_summary()):
            comment_parts.append(f"Criteria applied: {criteria_summary}")

        if comment_parts:
            full_comment = " ".join(comment_parts)
            full_comment = full_comment.replace("\n", " ")
            data["comment"] = full_comment

        if date_last_evaluated := self.value(SpecialEKeys.CURATION_DATE):
            data["dateLastEvaluated"] = self.convert_date(date_last_evaluated)
        elif date_last_reviewed := self.value(SpecialEKeys.CURATION_VERIFIED_DATE):
            data["dateLastEvaluated"] = self.convert_date(date_last_reviewed).with_more_messages(JsonMessages.warning(
                "No curation date provided, falling back to curation verified date."))
        return ValidatedJson(data)

    @property
    def germline_classification(self) -> ValidatedJson:
        data = {
            "germlineClassificationDescription": self.clinvar_value(SpecialEKeys.CLINICAL_SIGNIFICANCE).value(single=True)
        }
        return self.base_classification + ValidatedJson(data)

    @property
    def oncogenicity_classification(self) -> ValidatedJson:
        value = self.clinvar_value(SpecialEKeys.CLINICAL_SIGNIFICANCE).value(single=True)
        raw_value = value.json_data
        messages = JSON_MESSAGES_EMPTY
        if raw_value and raw_value not in {"Oncogenic", "Likely oncogenic", "Uncertain significance", "Likely benign", "Benign"}:
            messages += JsonMessages.error(f"Unsupported value for Oncogenic value - {value}")
        data = {"oncogenicityClassificationDescription": value}
        return self.base_classification + ValidatedJson(data, messages)

    @property
    def clinical_impact_classification(self) -> ValidatedJson:
        messages = JSON_MESSAGES_EMPTY
        value = self.clinvar_value(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).value(single=True)
        raw_value = value.json_data
        if raw_value in {"Tier I - Strong", "Tier II - Potential"}:
            messages += JsonMessages.error("ClinVar requires a detailed assertion type for Clinical Impact of Tier I or Tier II that we can't currently provide")
        data = {"clinicalImpactClassificationDescription": value}
        return self.base_classification + ValidatedJson(data, messages)

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
                    if join in (MultiCondition.CO_OCCURRING, MultiCondition.UNCERTAIN):
                        data["multipleConditionExplanation"] = join.clinvar_label
                    else:
                        messages += JsonMessages.error("ClinVar API only supports Co-occurring or Uncertain for multiple conditions")

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

        attempt_allele_origin_default = True
        allele_origin_value = ""
        allele_origin_messages = JSON_MESSAGES_EMPTY

        bucket: Optional[AlleleOriginBucket] = None
        if based_on := self.clinvar_export_record.classification_based_on:
            bucket = based_on.classification.allele_origin_bucket
            if bucket != AlleleOriginBucket.GERMLINE:
                attempt_allele_origin_default = False

        if allele_origin := self.clinvar_value(SpecialEKeys.ALLELE_ORIGIN).value(single=True, optional=True):
            allele_origin_value = allele_origin
        elif attempt_allele_origin_default:
            # default allele origin based on bucket if no value is provided (as there is another setting that can mean a blank all
            if bucket == AlleleOriginBucket.GERMLINE:
                allele_origin_value = "germline"
                allele_origin_messages += JsonMessages.info("Defaulting \"Allele origin\" to \"germline\" as no value provided")
            elif bucket == AlleleOriginBucket.SOMATIC:
                allele_origin_value = "somatic"
                allele_origin_messages += JsonMessages.info("Defaulting \"Allele origin\" to \"somatic\" as no value provided")

        data["alleleOrigin"] = ValidatedJson(allele_origin_value, allele_origin_messages)

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
