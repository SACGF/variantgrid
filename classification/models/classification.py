import copy
import re
import uuid
from collections import Counter, namedtuple
from dataclasses import dataclass
from datetime import datetime, timezone, timedelta, date
from enum import Enum, StrEnum
from functools import cached_property
from typing import Any, Dict, List, Union, Optional, Iterable, Callable, Mapping, TypedDict, Tuple, Set

import django.dispatch
from datetimeutc.fields import DateTimeUTCField
from dateutil.tz import gettz
from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db import models, transaction
from django.db.models import TextField
from django.db.models.deletion import CASCADE, PROTECT, SET_NULL
from django.db.models.expressions import RawSQL, OuterRef, Value, Subquery
from django.db.models.functions import LPad, Cast, Concat
from django.db.models.query import QuerySet
from django.db.models.query_utils import Q
from django.dispatch.dispatcher import receiver
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel
from guardian.shortcuts import assign_perm, get_objects_for_user
from annotation.models.models import AnnotationVersion, VariantAnnotationVersion, VariantAnnotation
from annotation.regexes import db_ref_regexes, DbRegexes
from classification.enums import ClinicalSignificance, SubmissionSource, ShareLevel, SpecialEKeys, \
    CRITERIA_NOT_MET, ValidationCode, CriteriaEvaluation, WithdrawReason, AlleleOriginBucket
from classification.models.classification_import_run import ClassificationImportRun
from classification.models.classification_patcher import patch_fuzzy_age
from classification.models.classification_utils import \
    ValidationMerger, ClassificationJsonParams, PatchMeta, ClassificationPatchResponse
from classification.models.classification_variant_info_models import ImportedAlleleInfo, ImportedAlleleInfoStatus
from classification.models.evidence_key import EvidenceKeyValueType, \
    EvidenceKey, EvidenceKeyMap, VCDataDict, WipeMode, VCDataCell, EvidenceKeyOverrides
from classification.models.evidence_mixin import EvidenceMixin, VCPatch
from classification.models.evidence_mixin_summary_cache import ClassificationSummaryCalculator, \
    ClassificationSummaryCacheDict
from classification.models.flag_types import classification_flag_types
from flags.models import Flag, FlagPermissionLevel, FlagStatus
from flags.models.models import FlagsMixin, FlagCollection, FlagTypeContext, \
    flag_collection_extra_info_signal, FlagInfos
from genes.hgvs import HGVSMatcher, CHGVS
from genes.models import Gene
from library.cache import clear_cached_property
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.guardian_utils import clear_permissions
from library.log_utils import report_exc_info, report_event
from library.preview_request import PreviewData, PreviewModelMixin, PreviewKeyValue
from library.utils import empty_to_none, nest_dict, cautious_attempt_html_to_text, \
    invalidate_cached_property, md5sum_str, get_timer
from ontology.models import OntologyTerm, OntologySnake, OntologyTermRelation
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Variant, Lab, Sample
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import AlleleSource, Allele, VariantAllele
from snpdb.user_settings_manager import UserSettingsManager

ChgvsKey = namedtuple('CHGVS', ['short', 'column', 'build'])

classification_validation_signal = django.dispatch.Signal()  # args: "classification", "patch_meta", "key_map"
classification_current_state_signal = django.dispatch.Signal()  # args: "user"
classification_post_publish_signal = django.dispatch.Signal()  # args: "classification", "previously_published", "previous_share_level", "newly_published", "user"
classification_withdraw_signal = django.dispatch.Signal()  # args: "classification", "user"
classification_variant_set_signal = django.dispatch.Signal()  # args: "classification", "variant"
classification_revalidate_signal = django.dispatch.Signal()  # args "classification"
# This signal is called if a classification is assigned to a variant (or liftover), or the classification changes
variants_classification_changed_signal = django.dispatch.Signal()  # args: "classification", "variant"

_key_to_regex = {
    'db_rs_id': DbRegexes.SNP,
    'uniprot_id': DbRegexes.UNIPROTKB,
    'clinvar_variation_id': DbRegexes.CLINVAR,
    'cosmic_id': DbRegexes.COSMIC,
    'gtr_id': DbRegexes.GTR
}


class VCBlobKeys(Enum):
    VALUE = "value"
    NOTE = "note"
    DB_REFS = "db_refs"
    EXPLAIN = "explain"


class ClassificationProcessError(Exception):
    """
    Use to report critical errors that an API user should be able to see
    e.g. referring to a non-existent lab
    """


class CreateNoClassificationForbidden(PermissionDenied):
    """ Used if settings prevent creating classifications not linked to a variant """


class ClassificationImport(models.Model):
    """
    A model to link together allele infos that will need variant matching.
    (Previously linked directly to Classifications, thus the name)
    """
    created = models.DateTimeField(auto_now_add=True)
    user = models.ForeignKey(User, on_delete=CASCADE)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)

    def get_variants_qs(self) -> QuerySet[Variant]:
        mv_ids = ImportedAlleleInfo.objects.filter(classification_import=self).values_list('matched_variant', flat=True)
        return Variant.objects.filter(pk__in=mv_ids)

    def __str__(self):
        return f"ClassificationImport ({self.genome_build})"


class ClassificationImportAlleleSource(AlleleSource):
    """ A model to indicate that variants need to be linked to an allele and lifted over to other builds """
    classification_import = models.ForeignKey(ClassificationImport, null=True, on_delete=CASCADE)

    def get_genome_build(self) -> GenomeBuild:
        return self.classification_import.genome_build

    def get_allele_qs(self) -> QuerySet[Allele]:
        variants_qs = self.classification_import.get_variants_qs()
        return Allele.objects.filter(variantallele__variant__in=variants_qs)

    def liftover_complete(self, genome_build: GenomeBuild):
        allele_qs = self.get_allele_qs()
        # Populate ClinGen for newly created variants as well as it's possible that old one failed but new could work
        build_variants = Variant.objects.filter(variantallele__genome_build=genome_build,
                                                variantallele__allele__in=allele_qs)
        populate_clingen_alleles_for_variants(genome_build, build_variants)

        ImportedAlleleInfo.relink_variants(self.classification_import)

        report_event('Completed import liftover',
                     extra_data={'liftover_id': self.pk, 'allele_count': allele_qs.count()})


class AllClassificationsAlleleSource(TimeStampedModel, AlleleSource):
    """ Used to reload all Classifications (for upgrades etc.) """

    script = models.TextField()  # Set if run from script, can check hasn't been re-run
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    git_hash = models.TextField()

    def get_genome_build(self) -> GenomeBuild:
        return self.genome_build

    def get_variants_qs(self) -> QuerySet[Variant]:
        # Note: This deliberately only gets classifications where the submitting variant was against this genome build
        # ie we don't use Classification.get_variant_q_from_classification_qs() to get liftovers
        contigs_q = Variant.get_contigs_q(self.genome_build)
        return Variant.objects.filter(contigs_q, importedalleleinfo__isnull=False)

    def liftover_complete(self, genome_build: GenomeBuild):
        ImportedAlleleInfo.relink_variants()


@receiver(flag_collection_extra_info_signal, sender=FlagCollection)
def get_extra_info(flag_infos: FlagInfos, user: User, **kwargs) -> None:  # pylint: disable=unused-argument
    """
    Allows us to provide extra info for FlagCollections attached to Classification
    e.g. linking to the appropriate allele page, discordance report etc.
    :param flag_infos: Information on the flag collections being displayed to the user.
    Populates this with the extra info
    :param user: The current user
    :param kwargs: Required by @receiver
    """
    from classification.models.discordance_models import DiscordanceReportClassification
    from classification.enums.discordance_enums import DiscordanceReportResolution

    vcs = Classification.objects.filter(flag_collection__in=flag_infos.ids).select_related('lab')
    drcs = DiscordanceReportClassification.objects.filter(classification_original__classification__in=vcs,
                                                          report__resolution=DiscordanceReportResolution.ONGOING)
    drcs_dict = {}
    for drc in drcs.values_list('classification_original__classification', 'report'):
        drcs_dict[drc[0]] = drc[1]

    for vc in vcs:
        """
            # If we are to use sub-flags in their current incarnation, here's how we'd include
            # one if zygosity or condition are missing
            flags_for_vc = flag_infos.flags_for_collection[vc.flag_collection_id]
            unsubmitted_flag = find_unsubmitted(flags_for_vc)

            if unsubmitted_flag:
                if vc.get(SpecialEKeys.ZYGOSITY) is None:
                    flag_infos.record_sub_flag(unsubmitted_flag, 'Missing Zygosity', 'Z')
                if vc.get(SpecialEKeys.CONDITION) is None:
                    flag_infos.record_sub_flag(unsubmitted_flag, 'Missing Condition', 'C')
        """

        context = {
            'label': vc.friendly_label,
            'vc_id': vc.id
        }
        if vc.variant_id:
            context['variant'] = vc.variant_id
        if vc.clinical_context_id:
            context['clinical_context'] = vc.clinical_context_id
        if user and vc.can_write(user):
            context['can_write'] = True

        dc_report_id = drcs_dict.get(vc.id)
        if dc_report_id:
            context['discordance_report'] = dc_report_id

        flag_infos.set_extra_info(vc.flag_collection_id, context, source_object=vc)


class ConditionResolvedTermDict(TypedDict):
    term_id: str
    name: str


class ConditionResolvedDict(TypedDict, total=False):
    """
    Structure of data used to cached resolved condition text again a classification
    """
    display_text: str  # plain text to show to users if not in a position to render links
    sort_text: str  # lower case representation of description
    resolved_terms: list[ConditionResolvedTermDict]
    plain_text_terms: list[str]  # A full list of unresolved plain text conditions
    resolved_join: str


@dataclass(frozen=True)
class ConditionResolved:
    terms: List[OntologyTerm]
    plain_text_terms: List[str] = None
    join: Optional['MultiCondition'] = None
    plain_text: Optional[str] = None  # fallback, not populated in all contexts

    def __hash__(self):
        hash_total = 0
        if terms := self.terms:
            for t in terms:
                hash_total += hash(t)
        hash_total += hash(self.join, )
        hash_total += hash(self.plain_text)
        return hash_total

    @property
    def summary(self) -> str:
        text = ", ".join([term.id for term in self.terms])
        if join := self.join:
            try:
                join = join.value
            except Exception:
                pass
            text = f"{text} {join}"
        return text

    @staticmethod
    def from_dict(condition_dict: ConditionResolvedDict) -> 'ConditionResolved':
        terms = [OntologyTerm.get_or_stub_cached(term.get("term_id")) for term in condition_dict.get("resolved_terms")]
        join = None
        if len(terms) > 1:
            from classification.models import MultiCondition
            join = MultiCondition(condition_dict.get("resolved_join"))

        terms.sort()
        return ConditionResolved(
            terms=terms,
            plain_text_terms=condition_dict.get("plain_text_terms"),
            join=join
        )

    @property
    def is_multi_condition(self) -> bool:
        return len(self.terms) > 1

    @property
    def single_term(self) -> Optional[OntologyTerm]:
        """
        Only call if not is_multi_condition
        """
        return self.terms[0] if len(self.terms) == 1 else None

    @cached_property
    def mondo_term(self) -> Optional[OntologyTerm]:
        if term := self.single_term:
            return OntologyTermRelation.as_mondo(term)
        else:
            return None

    def as_mondo_if_possible(self) -> 'ConditionResolved':
        if mondo_term := self.mondo_term:
            return ConditionResolved(terms=[mondo_term])
        else:
            return self

    def is_same_or_more_specific(self, other: 'ConditionGroup') -> bool:
        if self.is_multi_condition or other.is_multi_condition:
            # when looking at multiple conditions, do not attempt merging unless we're the exact same
            return self.terms == other.terms and self.join == other.join
        elif self.single_term == other.single_term:
            return True
        else:
            if other_mondo := other.mondo_term:
                if self_mondo := self.mondo_term:
                    descendant_relationships = OntologySnake.check_if_ancestor(descendant=self_mondo,
                                                                               ancestor=other_mondo)
                    return bool(descendant_relationships)

            # terms cant be converted to MONDO and not exact match, just return False
            return False

    @staticmethod
    def more_general_term_if_related(resolved_1: 'ConditionResolved', resolved_2: 'ConditionResolved') -> Optional[
        'ConditionResolved']:
        more_general: Optional[ConditionResolved] = None
        if resolved_1.is_same_or_more_specific(resolved_2):
            more_general = resolved_2
        elif resolved_2.is_same_or_more_specific(resolved_1):
            more_general = resolved_1

        if more_general:
            # if presented with different types, and we can switch over to MONDO, do so
            if not more_general.is_multi_condition and \
                    resolved_1.single_term.ontology_service != resolved_2.single_term.ontology_service:
                if mondo_term := more_general.mondo_term:
                    more_general = ConditionResolved(terms=[mondo_term])
            return more_general

        return None

    @staticmethod
    def term_to_dict(term: OntologyTerm) -> ConditionResolvedTermDict:
        term_dict: ConditionResolvedTermDict = {
            "term_id": term.id,
            "name": term.name
        }
        return term_dict

    @property
    def as_plain_text(self) -> str:
        if text := self.plain_text:
            return text
        else:
            return self.to_json()['display_text']

    def to_json(self, include_join: bool = True) -> ConditionResolvedDict:
        jsoned: ConditionResolvedDict
        if self.terms:
            from classification.models import MultiCondition

            def format_term(term: OntologyTerm) -> str:
                if name := term.name:
                    return f"{term.id} {name}"
                return term.id

            terms = self.terms
            text = ", ".join([format_term(term) for term in terms])
            if self.plain_text_terms:
                text += ",".join(self.plain_text_terms)

            sort_text = ", ".join([term.name for term in terms]).lower()
            join: Optional[MultiCondition] = None
            if len(terms) > 1 and include_join:
                join = self.join or MultiCondition.NOT_DECIDED
                text = f"{text}; {join.label}"

            resolved_term_dicts: List[ConditionResolvedTermDict] = [ConditionResolved.term_to_dict(term) for term in
                                                                    self.terms]
            jsoned: ConditionResolvedDict = {
                "resolved_terms": resolved_term_dicts,
                "resolved_join": join,
                "plain_text_terms": self.plain_text_terms,
                "display_text": text,
                "sort_text": sort_text
            }
            return jsoned
        else:
            jsoned: ConditionResolvedDict = {
                "plain_text_terms": self.plain_text_terms,
                "display_text": ", ".join(pt.lower() for pt in self.plain_text) if self.plain_text else None,
                "sort_text": ", ".join(pt.lower() for pt in self.plain_text) if self.plain_text else None
            }
        return jsoned

    @cached_property
    def join_text(self) -> Optional[str]:
        if join := self.join:
            try:
                from classification.models import MultiCondition
                return MultiCondition(join).label
            except:
                pass
        return None

    def __eq__(self, other: 'ConditionResolved') -> bool:
        if (s_terms := self.terms) and (o_terms := other.terms):
            return s_terms == o_terms and self.join == other.join
        elif self.terms or other.terms:
            return False
        else:
            return self.plain_text == other.plain_text

    def __lt__(self, other: 'ConditionResolved') -> bool:
        self_terms = self.terms or []
        other_terms = other.terms or []
        if len(self_terms) != len(other_terms):
            return len(self_terms) < len(other_terms)
        if len(self_terms) >= 1:
            return self_terms[0] < other_terms[0]
        return (self.plain_text or '') < (other.plain_text or '')


class ClassificationOutstandingIssues:

    def __init__(self, classification: 'Classification'):
        self.classification = classification
        self.flags = []
        self.issues = []

    @property
    def pk(self):
        return self.classification.pk

    def add_flag(self, flag: str):
        self.flags.append(flag)

    def add_issue(self, issue: str):
        self.issues.append(issue)

    def __str__(self):
        return f"({self.classification.friendly_label}) {', '.join(self.issues)} {', '.join(self.flags)}"


class Classification(GuardianPermissionsMixin, FlagsMixin, EvidenceMixin, TimeStampedModel, PreviewModelMixin):
    """
    A Variant Classification, belongs to a lab and user. Keeps a full history using ClassificationModification
    The data is free form basked on EvidenceKey (rather than one column per possible field)

    Based on ACMG criteria and evidence - @see https://www.nature.com/articles/gim201530
    """

    # TODO - remove variant and allele in favour of having that accessed via  variant_info
    variant = models.ForeignKey(Variant, null=True, on_delete=PROTECT)  # Null as might not match this
    """ Deprecated -  """
    allele = models.ForeignKey(Allele, null=True, on_delete=PROTECT)

    allele_info = models.ForeignKey(ImportedAlleleInfo, null=True, blank=True, on_delete=SET_NULL)
    """ Keeps links to common builds (37, 38) for quick access to c.hgvs, transcript etc. Is shared between 
        classifications with same import data """

    @property
    def allele_object(self) -> Allele:
        """ The new preferred way to reference the allele, so we can eventually remove allele from the
            classification object """
        try:
            return self.allele_info.allele
        except AttributeError:
            return self.allele

    sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL)
    """ The sample (if any) this classification is linked to """

    # classification_import = models.ForeignKey(ClassificationImport, null=True, on_delete=CASCADE)
    """
    provide a value here to later match this record to a variant
    classification_import has been removed and moved to allele_info
    """

    user = models.ForeignKey(User, on_delete=PROTECT)
    """ The owner of the classification, is somewhat redundant to the evidence_key classified_by, not heavily used """

    lab = models.ForeignKey(Lab, on_delete=CASCADE)
    """ The lab the classification belongs to, determines what permissions the classification is given """

    share_level = models.CharField(max_length=16, choices=ShareLevel.choices(), null=False, default=ShareLevel.LAB.key)
    """ The current share level of the classification, combined with lab determines the permissions """

    annotation_version = models.ForeignKey(AnnotationVersion, null=True, blank=True, on_delete=SET_NULL)
    """ If created from a variant and auto-populated, with which version of annotations. If null was
        created via import """

    clinical_context = models.ForeignKey('ClinicalContext', null=True, blank=True, on_delete=SET_NULL)
    """ After being matched to a variant, this will be set to the default clinical_context for the allele
    but it can be changed to be another clinical_context (for the same allele) """

    lab_record_id = models.TextField(blank=True, null=True)
    """ Should be unique together with lab """

    evidence = models.JSONField(null=False, blank=True, default=dict)
    """ The latest evidence (should always match the content of the latest ClassificationModification.evidence) """

    withdrawn = models.BooleanField(default=False)
    """ Soft delete, if withdrawn classification should not appear in most places """

    withdraw_reason = models.CharField(max_length=50, choices=WithdrawReason.choices, null=True, blank=True)
    """ Reason for withdrawing the classification """

    clinical_significance = models.CharField(max_length=1, choices=ClinicalSignificance.CHOICES, null=True, blank=True)
    """ Used as an optimisation for queries, is relatively out of date now """

    allele_origin_bucket = models.CharField(max_length=1, choices=AlleleOriginBucket.choices,
                                            default=AlleleOriginBucket.GERMLINE)
    """ Used to cache if we consider this classification germline or somatic """

    condition_resolution = models.JSONField(null=True, blank=True)  # of type ConditionProcessedDict

    summary = models.JSONField(null=False, blank=True, default=dict)  # useful for overall classification details
    """ Will be a ClassificationSummaryCacheDict """

    @property
    def summary_typed(self) -> ClassificationSummaryCacheDict:
        return self.summary

    last_source_id = models.TextField(blank=True, null=True)
    last_import_run = models.ForeignKey(ClassificationImportRun, null=True, blank=True, on_delete=SET_NULL)

    class Meta:
        unique_together = ('lab', 'lab_record_id')
        indexes = [
            models.Index(fields=["share_level"]),
            models.Index(fields=["withdrawn"]),
            models.Index(fields=["allele_origin_bucket"]),
            models.Index(models.F("summary__pathogenicity__sort"), name="summary__p_sort_idx"),
            models.Index(models.F("summary__somatic__sort"), name="summary__s_sort_idx"),
            models.Index(models.F("summary__date__value"), name="summary__d_sort_idx")
        ]

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-clipboard"

    @property
    def preview(self) -> PreviewData:
        title: Optional[str] = None
        extras = []
        if last_published := self.last_published_version:
            allele_origin_bucket = self.allele_origin_bucket_obj
            extras.append(
                PreviewKeyValue(
                    key="Allele Origin Grouping",
                    value=allele_origin_bucket.label
                )
            )
            classification_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            extras.append(PreviewKeyValue(
                key=classification_key.pretty_label,
                value=classification_key.pretty_value_from(last_published, empty_value="No Data")
            ))
            if allele_origin_bucket != AlleleOriginBucket.GERMLINE:
                somatic_clin_sig = EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)
                extras.append(PreviewKeyValue(
                    key=somatic_clin_sig.pretty_label,
                    value=somatic_clin_sig.pretty_value_from(last_published, empty_value="No Data")
                ))
        else:
            title = "Never been published"
        return self.preview_with(
            identifier=self.friendly_label,
            title=title,
            summary_extra=extras
        )

    @staticmethod
    def is_supported_transcript(transcript_or_hgvs: str):
        return ImportedAlleleInfo.is_supported_transcript(transcript_or_hgvs)

    @property
    def include_based_on_allele_info(self) -> bool:
        if allele_info := self.allele_info:
            if validation := allele_info.latest_validation:
                return validation.include
        return False

    @property
    def chgvs_grch37(self) -> Optional[str]:
        try:
            return self.allele_info.grch37.c_hgvs
        except AttributeError:
            return None

    @property
    def chgvs_grch38(self) -> Optional[str]:
        try:
            return self.allele_info.grch38.c_hgvs
        except AttributeError:
            return None

    @property
    def chgvs_grch37_compat(self) -> Optional[str]:
        try:
            return self.allele_info.grch37.c_hgvs_compat
        except AttributeError:
            return None

    @property
    def chgvs_grch38_compat(self) -> Optional[str]:
        try:
            return self.allele_info.grch38.c_hgvs_compat
        except AttributeError:
            return None

    @property
    def metrics_logging_key(self) -> Tuple[str, Any]:
        return "classification_id", self.pk

    @property
    def clinical_grouping_name(self) -> str:
        if self.share_level not in ShareLevel.DISCORDANT_LEVEL_KEYS:
            return 'not-shared'
        if cc := self.clinical_context:
            return cc.allele_origin_bucket + " " + cc.name
        return 'not-matched'

    @property
    def condition_text_record(self) -> 'ConditionText':
        try:
            if ctm := self.conditiontextmatch:
                return ctm.condition_text
        except Classification.conditiontextmatch.RelatedObjectDoesNotExist:
            pass
        return None

    @property
    def condition_resolution_dict(self) -> Optional[ConditionResolvedDict]:
        return self.condition_resolution

    @cached_property
    def condition_resolution_obj(self) -> Optional[ConditionResolved]:
        if cr_dict := self.condition_resolution_dict:
            return ConditionResolved.from_dict(cr_dict)
        return None

    def refresh_condition_resolution_details(self) -> bool:
        """
        If the condition sort text / display text has gotten out of date (as the term has updated itself) this will
        update the record
        :return: True if a change was detected and saved, False if no change was detected
        """
        if cond_obj := self.condition_resolution_obj:
            cond_json = cond_obj.to_json()
            if self.condition_resolution != cond_json:
                self.condition_resolution = cond_json
                invalidate_cached_property(self, 'condition_resolution_obj')
                self.save()
                return True
        return False

    @staticmethod
    def can_create_via_web_form(user: User):
        can_create = settings.CLASSIFICATION_WEB_FORM_CREATE
        return can_create and (user.is_superuser or settings.CLASSIFICATION_WEB_FORM_CREATE_BY_NON_ADMIN)

    @staticmethod
    def dashboard_report_new_classifications(since) -> int:
        return Classification.objects.filter(created__gte=since).count()

    @staticmethod
    def dashboard_total_shared_classifications() -> int:
        return Classification.objects.filter(lab__external=False, share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS,
                                             withdrawn=False).exclude(lab__name__icontains='legacy').count()

    @staticmethod
    def dashboard_total_unshared_classifications() -> int:
        qs = Classification.objects.filter(lab__external=False, withdrawn=False).exclude(lab__name__icontains='legacy')
        return qs.exclude(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).count()

    @staticmethod
    def dashboard_report_classifications_of_interest(since) -> List[ClassificationOutstandingIssues]:
        min_age = datetime.utcnow().replace(tzinfo=timezone.utc) - timedelta(
            minutes=2)  # give records 2 minutes to matching properly before reporting

        time_range_q = Q(created__gte=since) & Q(created__lte=min_age)

        # want to find new tags that are still open
        flag_collections = Flag.objects.filter(time_range_q, resolution__status=FlagStatus.OPEN)
        flag_collections = flag_collections.order_by('collection__id').values_list('collection__id', flat=True)
        flag_q = Q(flag_collection_id__in=flag_collections.distinct())
        missing_chgvs_q = (Q(allele_info__grch37__c_hgvs__isnull=True) | Q(allele_info__grch38__c_hgvs__isnull=True))
        coi_qs = Classification.objects.filter(flag_q | (time_range_q & missing_chgvs_q))
        coi_qs = coi_qs.order_by('-pk').select_related('lab', 'flag_collection')

        summaries: List[ClassificationOutstandingIssues] = []
        c: Classification
        for c in coi_qs:
            coi = ClassificationOutstandingIssues(c)
            this_flags = Flag.objects.filter(time_range_q, resolution__status=FlagStatus.OPEN,
                                             collection=c.flag_collection)

            variant_matching = False
            for flag_type in this_flags.order_by('flag_type').values_list('flag_type', flat=True):
                if flag_type == 'classification_matching_variant':
                    variant_matching = True
                coi.add_flag(flag_type)
            if not variant_matching:
                if not c.chgvs_grch37:
                    coi.add_issue("No cached 37 representation")
                if not c.chgvs_grch38:
                    coi.add_issue("No cached 38 representation")

            summaries.append(coi)

        return summaries

    @classmethod
    def order_by_evidence(cls, key_id: str):
        return RawSQL('cast(evidence->>%s as jsonb)->>%s', (key_id, 'value'))

    @property
    def variant_coordinate(self):
        """
        Used for the admin screen, so we can show variant coordinate in listing
        """
        return self.get(SpecialEKeys.VARIANT_COORDINATE)

    @property
    def imported_genome_build(self):
        return self.get(SpecialEKeys.GENOME_BUILD)

    @property
    def imported_c_hgvs(self) -> str:
        if c_hgvs := self.get(SpecialEKeys.C_HGVS):
            # remove any white space inside the c.HGVS
            c_hgvs = re.sub(r'\s+', '', c_hgvs)
            return c_hgvs

    @property
    def imported_g_hgvs(self):
        return self.get(SpecialEKeys.G_HGVS)

    def calc_allele_origin_bucket(self) -> AlleleOriginBucket:
        return AlleleOriginBucket.bucket_for_allele_origin(self.get(SpecialEKeys.ALLELE_ORIGIN))

    @property
    def allele_origin_bucket_obj(self) -> AlleleOriginBucket:
        return AlleleOriginBucket(self.allele_origin_bucket)

    def flag_type_context(self):
        return FlagTypeContext.objects.get(pk='classification')

    def flag_user_permission(self, user: User) -> FlagPermissionLevel:
        if self.can_write(user):
            return FlagPermissionLevel.OWNER

        # view permission is on the modification (not the variant classification)
        # but maybe it should be on both ?
        lp = self.last_published_version
        if lp and lp.can_view(user):
            return FlagPermissionLevel.USERS

        return FlagPermissionLevel.NO_PERM

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.requires_auto_population = False

    # Here to make this interchangeable with a ClassificationRef
    @property
    def record(self):
        return self

    @property
    def _evidence(self):
        return self.evidence

    @property
    def id_str(self):
        return str(self.id)

    @property
    def cr_lab_id(self):
        if settings.CLASSIFICATION_ID_OVERRIDE_PREFIX:
            return f"CR_{self.id}"
        return self.lab_record_id

    @property
    def friendly_label(self):
        return self.lab.name + ' / ' + self.cr_lab_id

    @staticmethod
    def to_date(date_str: str) -> Optional[datetime]:
        if date_str:
            return datetime.strptime(date_str, '%Y-%m-%d')
        return None

    @staticmethod
    def to_date_str(datetime_value: Union[date, datetime]) -> str:
        return datetime_value.strftime('%Y-%m-%d')

    def set_withdrawn(self, user: User, withdraw: bool = False, reason: str = 'OTHER') -> bool:
        if not self.id and withdraw:
            raise ValueError('Cannot withdrawn new classification record - use delete instead')

        if self.withdrawn == withdraw:
            return False  # no change

        self.withdrawn = withdraw
        self.withdraw_reason = reason if withdraw else None
        self.save()
        if withdraw:
            self.flag_collection_safe.get_or_create_open_flag_of_type(
                flag_type=classification_flag_types.classification_withdrawn,
                user=user,
                permission_check=False,
                reopen=True,
                comment=self.get_withdraw_reason_display()
            )
        else:
            self.flag_collection_safe.close_open_flags_of_type(
                flag_type=classification_flag_types.classification_withdrawn,
                user=user
            )

        classification_withdraw_signal.send(Classification, classification=self)
        return True

    def update_cached_c_hgvs(self):
        """
        :return: Returns length of the c.hgvs if successfully updated caches
        """
        self.update_allele()
        # self.update_allele_info_from_classification()
        self.allele_info.refresh_and_save(force_update=True)

    def update_allele(self):
        # Updates the allele based off the variant
        # Warning, does not call save()
        allele: Optional[Allele] = None
        if variant := self.variant:
            allele = variant.allele
        self.allele = allele

    def ensure_allele_info(self) -> Optional[ImportedAlleleInfo]:
        return self.ensure_allele_info_with_created()[0]

    def ensure_allele_info_with_created(self, force_allele_info_update_check: bool = False) -> Tuple[
        Optional[ImportedAlleleInfo], bool]:
        created = False
        if not self.allele_info or force_allele_info_update_check:
            try:
                genome_build_patch_version = self.get_genome_build_patch_version()
            except Exception:
                raise
                # no allele info if we can't derive a genome build
                # return None, False

            fields = {"imported_genome_build_patch_version": genome_build_patch_version}
            if c_hgvs := self.imported_c_hgvs:
                fields["imported_c_hgvs"] = c_hgvs
                fields["imported_md5_hash"] = md5sum_str(c_hgvs)
            elif g_hgvs := self.imported_g_hgvs:
                fields["imported_g_hgvs"] = g_hgvs
                fields["imported_md5_hash"] = md5sum_str(g_hgvs)
                fields["imported_transcript"] = self.transcript
            else:
                raise ValueError("Classification does not have imported_c_hgvs or imported_g_hgvs")
                # need either c.hgvs, g.hgvs (in the future, re-support variant_coordinate)
                # return None, False

            allele_info = ImportedAlleleInfo.get_or_create(**fields)
            if self.allele_info == allele_info:
                created = False
            else:
                self.allele_info = allele_info
                created = True
        return self.allele_info, created

    def update_allele_info_from_classification(self, force_update: bool = False) -> bool:
        """
        DEPRECATED
        Only to be called during migration from classification taking care of the matching to allele_info having the data
        """
        if self.allele_info:
            if not force_update and self.allele_info.status == ImportedAlleleInfoStatus.MATCHED_ALL_BUILDS:
                return False

        if allele_info := self.ensure_allele_info():
            # only need update_variant_coordinate for systems that were migrated when half of the AlleleInfo was done
            allele_info.update_variant_coordinate()
            allele_info.update_status()
            allele_info.apply_validation()
            allele_info.save()
            if not force_update and self.allele_info.status == ImportedAlleleInfoStatus.MATCHED_ALL_BUILDS:
                return False
            if matched_variant := self.variant:
                allele_info.set_variant_and_save(matched_variant, force_update=force_update)
            return True
        else:
            return False

    def attempt_set_variant_info_from_pre_existing_imported_allele_info(self) -> ImportedAlleleInfo:
        """
        Link to ImportedAlleleInfo (if we haven't already), then update classification with any existing ImportedAlleleInfo
        :return: True if there's nothing more to do, False if this may still require matching
        """
        allele_info = self.ensure_allele_info()
        self.apply_allele_info_to_classification()
        return allele_info

    def apply_allele_info_to_classification(self):
        if not self.allele_info:
            raise ValueError("Can't apply_allele_info when no allele_info has been set")

        allele_or_variant_changed = False
        variant: Optional[Variant] = None
        if vi := self.allele_info.variant_info_for_imported_genome_build:
            if viv := vi.variant:
                allele_or_variant_changed = self.variant != viv
                variant = viv
        allele_or_variant_changed = allele_or_variant_changed or self.allele_info

        if allele_or_variant_changed:
            self.allele = self.allele_info.allele
            self.variant = variant
            self.save()

    @property
    def share_level_enum(self) -> ShareLevel:
        return ShareLevel(self.share_level)

    @share_level_enum.setter
    def share_level_enum(self, value: ShareLevel):
        self.share_level = value.value

    @classmethod
    def filter_for_user(cls, user: User, queryset: Optional[QuerySet] = None, **kwargs) -> QuerySet:
        """ Classification only has write permission, View is based on a version in a point at time
            see ClassificationModification's read permission """
        klass = queryset if queryset is not None else cls

        if user and user.is_authenticated:
            queryset = get_objects_for_user(user, cls.get_write_perm(), klass=klass, accept_global_perms=True)
        else:
            queryset = cls.objects.none()

        return queryset

    @staticmethod
    def summarize_evidence_weights(evidence: dict, ekeys: Optional[EvidenceKeyMap] = None) -> str:
        if ekeys is None:
            ekeys = EvidenceKeyMap.instance()
        evidence_strings = []
        weights = {}
        for ekey in ekeys.all_keys:
            if ekey.value_type == EvidenceKeyValueType.CRITERIA:
                if ekey.key in evidence and evidence[ekey.key]:
                    value = evidence[ekey.key]
                    if isinstance(value, Mapping):
                        value = value.get('value')
                    if value == CRITERIA_NOT_MET:
                        # don't report on not-mets
                        continue
                    weights[value] = weights.get(value, 0) + 1

        for strength in CriteriaEvaluation.ALL_STRENGTHS:
            count = weights.get(strength, 0)
            if count > 0:
                # UNSPECIFIED STRENGTH HANDLING
                strength_str = strength
                if strength_str == CriteriaEvaluation.BENIGN_UNSPECIFIED:
                    strength_str = "B(Unspecified)"
                elif strength_str == CriteriaEvaluation.PATHOGENIC_UNSPECIFIED:
                    strength_str = "P(Unspecified)"
                evidence_strings.append(str(count) + 'x' + strength_str)

        return ', '.join(evidence_strings)

    @staticmethod
    def create_with_response(user: User,
                             lab: Lab,
                             lab_record_id: Optional[str] = None,
                             data: Optional[Dict[str, Any]] = None,
                             save: bool = True,
                             source: SubmissionSource = SubmissionSource.VARIANT_GRID,
                             make_fields_immutable=False,
                             populate_with_defaults=False,
                             **kwargs) -> Tuple['Classification', ClassificationPatchResponse]:
        """
            :param user: The user creating this Classification
            :param lab: The lab the record will be created under
            :param lab_record_id: The lab specific id for this record, will be auto-generated if not provided
            :param data: Initial data to populate with
            :param save: Should be True unless we're in test mode
            :param source: see patch_value
            :param make_fields_immutable: see patch_value
            :param populate_with_defaults: populate data with defaults as defined by the EvidenceKeys overrides
        """
        if data is None:
            data = {}

        making_new_id = False
        if lab_record_id is None:
            making_new_id = True
            lab_record_id = str(uuid.uuid4())  # temporary id
        record = Classification(user=user,
                                lab=lab,
                                lab_record_id=lab_record_id,
                                **kwargs)

        if populate_with_defaults:
            for e_key in record.evidence_keys.all_keys:
                if e_key.default_value is not None and e_key.key not in data:
                    data[e_key.key] = e_key.default_value

        response = record.patch_value(data,
                                      user=user,
                                      clear_all_fields=True,
                                      source=source,
                                      save=save,
                                      make_patch_fields_immutable=make_fields_immutable,
                                      initial_data=True,
                                      # revalidate all - so validations that only occur when specific fields change will still fire
                                      # even if no values are provided for that field, handy if there's validation that requires
                                      # at least 1 of 2 fields to be present for example
                                      revalidate_all=True)
        if save:
            try:
                if making_new_id:
                    record.lab_record_id = 'vc' + str(record.id)
                    record.save()
            except:
                pass
                # in case of collision, leave the uuid

        return record, response

    @staticmethod
    def create(user: User,
               lab: Lab,
               lab_record_id: Optional[str] = None,
               data: Optional[Dict[str, Any]] = None,
               save: bool = True,
               source: SubmissionSource = SubmissionSource.VARIANT_GRID,
               make_fields_immutable=False,
               populate_with_defaults=False,
               **kwargs) -> 'Classification':
        """
        Deprecated use create_with_response instead as otherwise so
        """
        return Classification.create_with_response(
            user=user,
            lab=lab,
            lab_record_id=lab_record_id,
            data=data,
            save=save,
            source=source,
            make_fields_immutable=make_fields_immutable,
            populate_with_defaults=populate_with_defaults,
            **kwargs
        )[0]

    def save(self, *args, **kwargs):
        """ The post_save event is fired after save(), for other logic that depends on Classifications
            but doesn't work with the Classification object itself.
        """
        fix_permissions = kwargs.pop('fix_permissions', None) or self.id is None
        super().save(*args, **kwargs)

        # fix permissions removes previous permissions
        if fix_permissions:
            self.fix_permissions()

    def get_key_errors(self) -> Dict[str, Dict]:
        """
        Returns a dict of key to validation error (first error in the case
        of multiple errors for one key).
        Only entries with at least one validation error will appear in the dict
        """
        key_errors = {}
        for key, entry in self.evidence.items():
            validation = entry.get('validation', [])
            for valid in validation:
                if valid.get('severity') == 'error':
                    key_errors[key] = valid
                    break
        return key_errors

    @staticmethod
    def _export_value(normal_value_obj: Any, key_data: EvidenceKey, export_key: str = 'key') -> Any:
        if options := key_data.virtual_options:
            if normal_value_obj and isinstance(normal_value_obj, Mapping):
                value = normal_value_obj.get('value')
            else:
                value = normal_value_obj

            for option in options:
                if export_key in option:
                    export_value = option.get(export_key)
                    if option.get('key') == value:
                        return export_value
        return None

    @staticmethod
    def match_option(options, check_value) -> Optional[Dict[str, str]]:
        for option in options:
            option_value = option.get('key')

            if check_value is True:
                if option.get('default'):
                    return option
                continue
            if not isinstance(check_value, str):
                check_value = str(check_value)

            aliases = []
            if option_value:
                aliases.append(option_value)
            if option.get('label'):
                aliases.append(option.get('label'))
            # Originally had index as a value you could submit, but this starts to get dangerous for fields
            # that allow custom values, can't distinguish between index 3 and custom value 3
            # if option.get('index') is not None:
            #    aliases.append(str(option.get('index')))
            elif option_value:
                aliases.append(EvidenceKey.pretty_label_from_string(option_value))

            aliases.extend(option.get('aliases', []))

            match = next((alias for alias in aliases if alias.lower() == check_value.lower()), None)
            if match is not None:
                return option
            check_value = check_value.replace('-', '_').replace(' ', '_')
            match = next((alias for alias in aliases if alias.lower() == check_value.lower()), None)
            if match is not None:
                return option
        return None

    @staticmethod
    def process_option_values(cell: VCDataCell, values: List[Any]) -> Optional[List[str]]:
        e_key = cell.e_key
        options = e_key.virtual_options or []
        # Do a case-insensitive check for each value against the key and any aliases
        # if there's a match to any of those, normalise back to the key (with the case of the key)
        results: List[str] = []
        # remove duplicates
        values = list(set(values))
        for check_value in values:
            check_value = empty_to_none(check_value)
            if check_value is None:
                continue
            if check_value is True and e_key.default_crit_evaluation:
                results.append(e_key.default_crit_evaluation)
                continue
            if check_value is False and e_key.value_type == EvidenceKeyValueType.CRITERIA:
                results.append(CriteriaEvaluation.NOT_MET)
                continue

            check_value = str(check_value)
            matched_option = Classification.match_option(options, check_value)
            if matched_option:
                results.append(matched_option.get('key'))
                if settings.CLASSIFICATION_REQUIRE_OVERWRITE_NOTE and \
                        matched_option.get('override') and cell.note is None and cell.explain is None:
                    cell.add_validation(
                        code=ValidationCode.REQUIRES_NOTE,
                        severity='warning',
                        message='Override value has been selected, but no note provided'
                    )
            else:
                # add the value we couldn't match
                results.append(str(check_value))

                if not e_key.allow_custom_values:
                    cell.add_validation(
                        code=ValidationCode.INVALID_VALUE,
                        severity='error',
                        message='Illegal value (' + str(check_value) + ')'
                    )

        if results:
            results.sort()  # TODO, sort by presentation index
            return results
        return None

    @cached_property
    def evidence_keys(self) -> EvidenceKeyMap:
        # note that this should be invalidated if the config changes (which would happen if
        # assertion method is updated for example)
        return EvidenceKeyMap.instance().with_overrides(self.evidence_key_overrides)

    def process_entry(self, cell: VCDataCell, source: str):
        """
        Normalises the value and adds any validation
        """
        # parse all notes looking for references
        value = cell.value
        e_key = cell.e_key
        note = cell.note

        if isinstance(value, str):
            if '\u00a0' in value:
                value = value.replace('\u00a0', ' ')
        if note and '\u00a0' in note:
            note = note.replace('\u00a0', ' ')

        if self.lab.external:
            cell.validate = False

        if value and e_key.value_type in (EvidenceKeyValueType.FREE_ENTRY, EvidenceKeyValueType.TEXT_AREA):
            scan_value = str(value)
            if results := db_ref_regexes.search(scan_value, default_regex=_key_to_regex.get(e_key.key)):
                cell.db_refs = [r.to_json() for r in results]

        if note:
            if results := db_ref_regexes.search(note):
                cell.db_refs = (cell.db_refs or []) + [r.to_json() for r in results]

        if source == SubmissionSource.API and isinstance(value, str):
            # we often get escaped > sent to us as part of
            value = value.replace('&gt;', '>').replace('&lt;', '<')

        if value is None and e_key.mandatory and not e_key.exclude_namespace:
            if e_key.key == SpecialEKeys.CONDITION and self.condition_resolution_obj:
                # Gene condition is not mandatory if we have a condition resolution
                pass
            else:
                # only provide mandatory validation for internal labs
                # don't want to dirty up imported read only classifications with errors
                cell.add_validation(code=ValidationCode.MANDATORY, severity='error', message='Missing mandatory value')

        # if we got an array of values
        if isinstance(value, list):
            if not value:
                # if the array is empty, treat that as a null value
                value = None
            elif e_key.value_type != EvidenceKeyValueType.MULTISELECT and len(value) == 1:
                # if we're not a multi-select key, and we have an array with 1 value in it
                # use that 1 value
                value = value[0]

        if value is not None:
            # remove support for existing_variant_id
            # variant matching is best done through allele_info stuff now
            if e_key.key == 'existing_sample_id':
                try:
                    self.sample = Sample.objects.get(pk=int(value))
                    # self.requires_auto_population = True
                    return None  # clear out the value
                except:
                    cell.add_validation(code=ValidationCode.MATCHING_ERROR, severity='error',
                                        message="Couldn't resolve Sample " + str(value) + ")")

            # convert users (username or email) to username if possible
            if e_key.value_type == EvidenceKeyValueType.USER:
                value = str(value)
                user = User.objects.filter(Q(username=value) | Q(email=value)).first()
                if user:
                    cell.value = user.username
                else:
                    pass

            # normalise booleans if it's un-ambiguous what the value is meant to be
            # otherwise leave the value as is to be caught by validation
            elif e_key.value_type == EvidenceKeyValueType.BOOLEAN:
                if isinstance(value, list) and len(value) == 1:
                    value = value[0]
                if not isinstance(value, bool):
                    if value == 1:
                        value = True
                    elif value == 0:
                        value = False
                    elif isinstance(value, str):
                        value = value.upper()
                        if value in ('TRUE', 'Y', 'YES', 'T', '1'):
                            value = True
                        elif value in ('FALSE', 'N', 'NO', 'F', '0'):
                            value = False
                        else:
                            cell.add_validation(code=ValidationCode.INVALID_VALUE, severity='error',
                                                message="Invalid value for yes/no field (" + str(value) + ")")
                cell.value = value

            # handle selects with multiple values (comma sep string, array etc.)
            elif e_key.value_type in (EvidenceKeyValueType.MULTISELECT,
                                      EvidenceKeyValueType.SELECT,
                                      EvidenceKeyValueType.CRITERIA):
                parts: List[Any]
                if isinstance(value, str):
                    if e_key.value_type == EvidenceKeyValueType.MULTISELECT and '|_' in str(value):
                        parts = value.split('|_')
                    else:
                        parts = value.split(',')
                    parts = [p.strip() for p in parts if p.strip()]
                elif isinstance(value, list):
                    parts = value
                else:
                    parts = [value]
                parts = [part for part in parts if part is not None and part != ""]  # filter out any blanks

                # sets the value of cell to the valid options
                value = Classification.process_option_values(cell, parts)
                if value is not None:
                    if e_key.value_type != EvidenceKeyValueType.MULTISELECT:
                        value_len = len(value)
                        if value_len > 1:
                            cell.add_validation(code=ValidationCode.TOO_MANY_VALUES, severity='error',
                                                message="Can only have a single value for a select field")
                        elif value_len == 1:
                            value = value[0]
                            cell.value = value
                cell.value = value

            # base text tidy
            elif e_key.value_type in (EvidenceKeyValueType.FREE_ENTRY, EvidenceKeyValueType.TEXT_AREA):
                value = str(value)
                if e_key.value_type == EvidenceKeyValueType.FREE_ENTRY:
                    # strip out HTML from single row files to keep our data simple
                    # allow HTML in text areas
                    if not value.startswith("http"):  # makes beautiful soap angry thinking we want to go to the URL
                        value = cautious_attempt_html_to_text(value)

                cell.value = value

            # -- BASIC NUMBER NORMALISING, VALIDATING --
            elif e_key.value_type == EvidenceKeyValueType.UNIT:
                valid = False
                try:
                    value = float(value)
                    cell.value = value
                    valid = 0 <= value <= 1
                except:
                    pass
                if not valid:
                    cell.add_validation(code=ValidationCode.INVALID_UNIT, severity='warning',
                                        message='Value should be a number between 0 and 1')

            elif e_key.value_type == EvidenceKeyValueType.FLOAT:
                valid = False
                try:
                    value = float(value)
                    cell.value = value
                    valid = True
                except:
                    pass
                if not valid:
                    cell.add_validation(code=ValidationCode.INVALID_FLOAT, severity='warning',
                                        message='Value should be a number')

            elif e_key.value_type == EvidenceKeyValueType.INTEGER:
                valid = False
                try:
                    value = int(value)
                    cell.value = value
                    valid = True
                except:
                    pass
                if not valid:
                    cell.add_validation(code=ValidationCode.INVALID_INTEGER, severity='warning',
                                        message='Value should be a whole number')

            elif e_key.value_type == EvidenceKeyValueType.DATE:
                try:
                    Classification.to_date(value)
                except ValueError as ve:
                    message = "Invalid date (expected yyyy-mm-dd)"
                    cell.add_validation(code=ValidationCode.INVALID_DATE, severity='warning',
                                        message=message)

            # only complain about unknown keys if we have a value
            if e_key.is_dummy:
                cell.add_validation(code=ValidationCode.UNKNOWN_KEY, severity='warning',
                                    message='This key ' + str(e_key.key) + ' has not been registered')

        # SPECIAL HANDLING OF REFSEQ, basically looking for transcript version
        if e_key.key == SpecialEKeys.REFSEQ_TRANSCRIPT_ID:
            if isinstance(value, str):
                version_regex = re.compile('.*?[.][0-9]+')
                if not version_regex.match(value):
                    cell.add_validation(code=ValidationCode.UNKNOWN_TRANSCRIPT, severity='warning',
                                        message='Transcript should include version e.g. NM_001256799.1')

        elif e_key.key == SpecialEKeys.AGE:
            if settings.CLASSIFICATION_AUTOFUZZ_AGE:
                cell.value = patch_fuzzy_age(value)

        # ensure we have one non None value before returning the structure
        if not cell.has_data:
            cell.wipe(WipeMode.SET_NONE)

    def modification_at_timestamp(self, version_timestamp) -> 'ClassificationModification':
        dt = datetime.utcfromtimestamp(version_timestamp).replace(tzinfo=timezone.utc)
        vcm = ClassificationModification.objects.filter(classification=self, created__lte=dt).order_by(
            '-created').first()
        return vcm

    def has_outstanding_changes(self) -> bool:
        last_edit = self.last_edited_version
        last_publish = self.last_published_version
        submitted = False
        if last_edit and last_publish:
            submitted = last_edit.pk == last_publish.pk
        return not submitted

    def revalidate(self, user: User, migration_patch: Optional[VCPatch] = None):
        """
        Re-normalises all values and re-performs validation on mutable fields.
        Re-publishes if the latest version was published (don't want to make every record
        end up with unshared changes warning)
        """
        submitted = not self.has_outstanding_changes()

        use_evidence = self.evidence.copy()
        e_keys = self.evidence_keys
        for key_entry in e_keys.mandatory():
            if key_entry.key not in use_evidence:
                use_evidence[key_entry.key] = None

        raw_patch = {}
        for key, value in use_evidence.items():
            if value is None:
                raw_patch[key] = None
            else:
                raw_patch[key] = {
                    "value": value.get('value'),
                    "note": value.get('note'),
                    "explain": value.get('explain'),
                    "immutable": value.get('immutable')
                }
        if migration_patch:
            for key, value in migration_patch.items():
                raw_patch[key] = value

        self.patch_value(
            patch=raw_patch,
            source=SubmissionSource.VARIANT_GRID,
            user=user,
            save=True,
            revalidate_all=True
        )
        # migration fixes are so retrieve_established at this point, no need to re-run them
        # self.fix_migration_stuff(user)

        classification_revalidate_signal.send(sender=Classification, classification=self)

        if submitted:
            # have to re-retrieve last_edited_version as it's now the latest
            self.last_edited_version.publish(share_level=self.share_level_enum, user=user, vc=self)

    def fix_migration_stuff(self, user: User) -> None:
        """
        Ensure there's a published version, and that it's marked as such, that the variant matching flags are correct
        """
        modifications_for_self_qs = ClassificationModification.objects.filter(classification=self)
        if not modifications_for_self_qs.filter(is_last_published=True).count() == 1:
            # in case we have multiple last published records, this will  ensure we only have 1
            modifications_for_self_qs.filter(is_last_published=True).update(is_last_published=False)
            last_published = modifications_for_self_qs.filter(published=True).order_by('-created').first()
            if last_published:
                last_published.is_last_published = True
                last_published.save()
            else:
                first_modification = modifications_for_self_qs.order_by('created').first()
                if first_modification:
                    first_modification.publish(share_level=ShareLevel.LAB, user=user, vc=self)

        if not modifications_for_self_qs.filter(is_last_edited=True).count() == 1:
            modifications_for_self_qs.filter(is_last_edited=True).update(is_last_edited=False)
            last_edited = modifications_for_self_qs.order_by('-created').first()
            if last_edited:
                last_edited.is_last_edited = True
                last_edited.save()

        if settings.UNSHARED_FLAG_ENABLED and self.share_level_enum.key not in ShareLevel.DISCORDANT_LEVEL_KEYS:
            self.flag_collection_safe.get_or_create_open_flag_of_type(
                flag_type=classification_flag_types.unshared_flag
            )

        self.apply_allele_info_to_classification()
        self.fix_permissions()
        self.save()

    def patch_history(self, patcher: Callable[[PatchMeta], None]):
        """
        Patches current evidence, and historic modifications "delta" and "published_evidence" fields
        :param patcher: Code that modified PatchMeta in place through the patch methods
        """
        last_evidence = {}
        modification: ClassificationModification
        for modification in ClassificationModification.objects.filter(classification=self).order_by('created'):
            patch_meta = PatchMeta(patch=modification.delta, existing=last_evidence)
            last_evidence = copy.deepcopy(modification.evidence)
            patcher(patch_meta)

            modification.delta = dict(patch_meta.patch)
            if modification.published_evidence:
                patch_meta = PatchMeta(patch=modification.published_evidence, existing={})
                patcher(patch_meta)
                modification.published_evidence = dict(modification.published_evidence)
            modification.save(update_fields=['delta', 'published_evidence'])

        self.evidence = self.last_edited_version.evidence
        self.save(update_fields=['evidence'])

    @transaction.atomic()
    def patch_value(self,
                    patch: Dict[str, Any],
                    clear_all_fields: bool = False,
                    user: Optional[User] = None,
                    source: SubmissionSource = None,
                    leave_existing_values=False,
                    save=False,
                    make_patch_fields_immutable=False,
                    remove_api_immutable=False,
                    initial_data=False,
                    revalidate_all=False,
                    ignore_if_only_patching: Optional[Set[str]] = None,
                    patch_known_keys_only: Optional[bool] = None) -> ClassificationPatchResponse:
        """
            Creates a new ClassificationModification if the patch values are different to the current values
            Patching a value with the same value has no effect
            :param patch: The set of values we're patching
            :param clear_all_fields: If completely overwriting the classification
            :param user: The user patching these values (for security checks) - should always be provided
            :param source: See SubmissionSource, determines what's immutable
            :param leave_existing_values: Only update empty fields
            :param save: saves to the database (leave as False if test mode or going to do other changes)
            :param make_patch_fields_immutable: Make all fields updated in this patch immutable.
            :param remove_api_immutable: If True, immutability level (under variantgrid) is removed from all fields.
                                         Requires source: SubmissionSource.VariantGrid
            :param initial_data: if True, divides c.hgvs to
            :param revalidate_all: if True, runs validation over all fields we have, otherwise only the values being patched
            :param ignore_if_only_patching: if provided, if only these fields are different in the patch, don't both to patch anything
            :param patch_known_keys_only: If true, data that doesn't match up with evidence keys will be ignored (e.g. {"acmg:bp34": "BM"} will do nothing)
            :returns: A dict with "messages" (validation errors, warnings etc.) and "modified" (fields that actually changed value)
        """
        if patch_known_keys_only is None:
            patch_known_keys_only = settings.CLASSIFICATION_ALLOW_UNKNOWN_KEYS is False

        source = source or SubmissionSource.API

        patch_response = ClassificationPatchResponse()

        is_admin_patch = source.can_edit(SubmissionSource.VARIANT_GRID)
        key_dict: EvidenceKeyMap = self.evidence_keys

        use_evidence = VCDataDict(copy.deepcopy(self.evidence),
                                  evidence_keys=self.evidence_keys)  # deep copy so don't accidentally mutate the data

        patch = VCDataDict(data=EvidenceMixin.to_patch(patch),
                           evidence_keys=self.evidence_keys)  # the patch we're going to apply on-top of the evidence

        if patch_known_keys_only:
            unrecognised_keys: set[str] = set()
            for key in patch.data.keys():
                if self.evidence_keys.get(key).is_dummy:
                    unrecognised_keys.add(key)
            if unrecognised_keys:
                for key in unrecognised_keys:
                    del patch.data[key]
                sorted_bad_keys = ", ".join(sorted(unrecognised_keys))
                patch_response.append_warning("unrecognised_key",
                                              f"The following keys are not valid ({sorted_bad_keys})")

        # make sure gene symbol is uppercase
        # need to do it here because it might get used in c.hgvs
        gene_symbol_cell = patch[SpecialEKeys.GENE_SYMBOL]
        if gene_symbol_cell.provided:
            if gene_symbol := gene_symbol_cell.value:
                if isinstance(gene_symbol, str) and gene_symbol != gene_symbol.upper():
                    gene_symbol_cell.value = gene_symbol.upper()

        # remove all whitespace from c.HGVS
        if SpecialEKeys.C_HGVS in patch:
            c_parts_cell = patch[SpecialEKeys.C_HGVS]
            if c_hgvs := c_parts_cell.value:
                c_hgvs = re.sub(r'\s+', '', c_hgvs)
                c_parts_cell.value = c_hgvs

        if initial_data:
            # if c.hgvs contains other values (such as
            if SpecialEKeys.C_HGVS in patch:
                c_parts_cell = patch[SpecialEKeys.C_HGVS]
                c_hgvs = c_parts_cell.value
                if c_hgvs:
                    c_parts = CHGVS(full_c_hgvs=c_hgvs)
                    transcript = c_parts.transcript
                    gene_symbol = c_parts.gene

                    # upper case gene symbol if it's not already
                    if gene_symbol and gene_symbol != gene_symbol.upper():
                        gene_symbol = gene_symbol.upper()
                        c_parts_cell.value = c_parts.with_gene_symbol(gene_symbol).full_c_hgvs

                    if transcript:
                        transcript_key = SpecialEKeys.REFSEQ_TRANSCRIPT_ID
                        if transcript.startswith('ENST'):
                            transcript_key = SpecialEKeys.ENSEMBL_TRANSCRIPT_ID
                        elif transcript.startswith('LRG_'):
                            transcript_key = SpecialEKeys.LRG_ID

                        transcript_cell = patch[transcript_key]
                        if not transcript_cell.provided:
                            transcript_cell.value = transcript

                    if gene_symbol and not gene_symbol_cell.provided:
                        # if no gene symbol value provided, populate it from c.hgvs
                        gene_symbol_cell.value = gene_symbol

                    elif not gene_symbol and gene_symbol_cell.value:
                        # if gene symbol provided (but not in c.hgvs) inject it into it
                        c_parts = c_parts.with_gene_symbol(gene_symbol_cell.value)
                        c_parts_cell.value = c_parts.full_c_hgvs

        # if submitting via API treat null as {value:None, explain:None, notes:None} for known keys,
        # so we clear out any previous values but still retain immutability
        # for non-existent keys just leave as null (if it started as null)
        if source == SubmissionSource.API:
            submitted_keys = list(patch.keys())
            for key in submitted_keys:
                cell = patch[key]
                if cell.raw is None and not cell.e_key.is_dummy:
                    cell.wipe(WipeMode.ATTRIBUTES_TO_NONE)

        if clear_all_fields:
            for key in use_evidence.keys():
                if key not in patch:
                    patch[key].wipe(WipeMode.SET_NONE)
            # when making a new classification ensure we trigger mandatory fields
            for e_key in key_dict.mandatory():
                if e_key not in patch:
                    patch[e_key].wipe(WipeMode.SET_EMPTY)

        if remove_api_immutable:
            for key in use_evidence.keys():
                cell = use_evidence[key]
                if cell.immutability == SubmissionSource.API:
                    patch[key].immutability = None

        for cell in patch.cells():
            if cell.raw is not None:

                # if we submitted via the API, provided a new value and did not provide an explain
                # explicitly set explain to None so any previous explanation is removed
                if not cell.explain and source != SubmissionSource.FORM:
                    cell.explain = None

                if not is_admin_patch:
                    # strip out any attribute that can't be provided by the client
                    cell.strip_non_client_submission()

                # if providing a value that's marked as immutable on evidence keys
                # set the immutability level to Variant Grid
                if source in [SubmissionSource.VARIANT_GRID, SubmissionSource.API] and \
                        cell.value is not None and cell.e_key.immutable:
                    cell.immutability = SubmissionSource.VARIANT_GRID

                elif make_patch_fields_immutable and not cell.immutability and source != SubmissionSource.FORM:
                    cell.immutability = source

                # Don't set immutability to form as that's the lowest immutability level
                if cell.immutability == SubmissionSource.FORM:
                    cell.pop_immutability()

                # copy attributes from the existing into the patch (they'll later be removed by diff)
                # required for process_entry
                existing = use_evidence[cell.e_key]

                if not cell.has_change(existing) and not revalidate_all:
                    patch.ignore(cell.e_key.key)
                    continue

                if existing.raw is not None:
                    cell.merge_from(existing)

            self.process_entry(cell, source=source)

        # creating the record, if no owner was provided set the value to match current
        if not self.id:
            owner_cell = patch[SpecialEKeys.OWNER]
            if not owner_cell.provided and user:
                owner_cell.value = user.username

        debug_timer = get_timer()
        debug_timer.tick("Patch data normalisation")

        patch_meta = PatchMeta(patch=patch.data, existing=use_evidence.data, revalidate_all=revalidate_all)

        validations_received = classification_validation_signal.send(sender=Classification, classification=self,
                                                                     patch_meta=patch_meta, key_map=key_dict)
        multi_errors = ValidationMerger.union_send_responses(validations_received)
        multi_errors.apply(patch.data, use_evidence.data)

        # merge the patch with a lot of nulls if it's a POST
        # i.e. a "complete" patch

        # now only include the values in the patch that
        # are diff to existing values
        diffs_only_patch = {}
        diffs_only_patch_data = VCDataDict(data=diffs_only_patch, evidence_keys=self.evidence_keys)

        for cell in patch.cells():
            existing = use_evidence[cell.e_key]
            key = cell.e_key.key
            existing_immutability = existing.immutability or SubmissionSource.FORM
            patched = diffs_only_patch_data[cell.e_key]

            if cell == existing:
                pass
            elif existing.value is not None and leave_existing_values:
                pass
            elif existing.raw is None:
                patched.raw = cell.diff(None)
            elif cell.raw is None:
                if not source.can_edit(existing_immutability):
                    message = f"Cannot change immutable value for {key} from {existing.value} to blank"
                    patch_response.append_warning(key=key, code="immutable", message=message)
                    patched.wipe(WipeMode.POP)  # reject entire change if attempting to change immutable value
                else:
                    patched.wipe(WipeMode.SET_NONE)
            else:
                patched.raw = cell.diff(dest=existing, ignore_if_omitted={'immutable'})
                if ('value' in patched or 'explain' in patched) and not source.can_edit(existing_immutability):
                    msg = f"Cannot change immutable value or explain for {key} from {existing.value} to {patched.value}"
                    patch_response.append_warning(key=key, code="immutable", message=msg)
                    patched.wipe(WipeMode.POP)  # reject entire change if attempting to change immutable value
                else:

                    # check to see if immutable is in both the source and value
                    # if the original immutable level is higher or equal than the new one, remove it for redundancy
                    new_immutability = patched.immutability
                    if new_immutability is not None:
                        if not new_immutability.can_edit(existing_immutability):
                            patched.pop_immutability()

                    # if we ended up with no attributes to patch, don't patch with an empty dict
                    # could be two floats that are close enough to each other
                    if not patched.has_data:
                        patched.wipe(WipeMode.POP)

        apply_patch = bool(diffs_only_patch)
        if apply_patch and ignore_if_only_patching and ignore_if_only_patching.issuperset(diffs_only_patch.keys()):
            apply_patch = False

        if not apply_patch:
            # if this is the first submission, create a version regardless of if there's data or not
            apply_patch = self.id is None

        pending_modification: Optional[ClassificationModification]

        debug_timer.tick("Patch validation")

        # only make a modification if there's data to actually patch
        if apply_patch:
            patch_response.modified_keys = set(diffs_only_patch.keys())
            last_edited = self.last_edited_version

            # in some cases we can append the last version
            # we do this so editing on the web form over the period of 1 minute
            # makes 2 versions instead of one every 2 seconds
            if last_edited and last_edited.is_edit_appendable(user=user, source=source):
                pending_modification = last_edited
                # update last modification as it was recent enough
                delta = pending_modification.delta
                delta.update(diffs_only_patch)
                pending_modification.delta = delta

            else:
                if save and last_edited:
                    last_edited.is_last_edited = False
                    last_edited.save()

                # create a new modification
                pending_modification = ClassificationModification(
                    user=user,
                    share_level=ShareLevel.LAB.key,
                    source=source.value,
                    delta=diffs_only_patch,
                    is_last_edited=True
                )
            ###
            # ACTUALLY APPLY THE PATCH
            # update the actual evidence with normalised patch values
            ###
            Classification.patch_with(target=use_evidence.data, patch=diffs_only_patch, tidy_nones=True)

            self.evidence = use_evidence.data
            # classification is stored on the classification record and on the classification modification
            # (not sure if we actually use it for anything on the modification)

            # TODO this is in the wrong spot, it should be on publish
            clinical_significance_choice = self.calc_clinical_significance_choice()
            self.clinical_significance = clinical_significance_choice
            pending_modification.clinical_significance = clinical_significance_choice

            self.allele_origin_bucket = self.calc_allele_origin_bucket()

            message = 'Patched changed values for ' + ', '.join(sorted(diffs_only_patch.keys()))
            patch_response.append_warning(code="patched", message=message)

            # update the evidence keys in case that's changed
            # FIXME do some optimisation so we can see if it changed easily
            self._clear_evidence_key_cache()

            if save:
                self.save()
                patch_response.saved = True
                # have to save the modification after saving the record in case this
                # is the first record
                if pending_modification:
                    pending_modification.classification = self
                    pending_modification.save()
                    assign_perm(ClassificationModification.get_read_perm(), self.lab.group, pending_modification)
            debug_timer.tick("Saved")

        # end apply patch diff
        else:
            if diffs_only_patch:
                message = 'Only ' + ', '.join(sorted(diffs_only_patch.keys())) + ' changed, ignoring'
                patch_response.append_warning(code="no_patch", message=message)
            else:
                patch_response.append_warning(code="no_patch", message="No changes detected to patch")
                # don't save if we haven't changed any values

        if self.requires_auto_population:
            self.requires_auto_population = False
            from classification.autopopulate_evidence_keys.autopopulate_evidence_keys import \
                classification_auto_populate_fields
            genome_build = self.get_genome_build()
            patch_response += classification_auto_populate_fields(self, genome_build, save=save)

        return patch_response

    def publish_latest(self, user: User, share_level=None) -> bool:
        if not share_level:
            share_level = self.share_level_enum
        latest_edited = self.last_edited_version
        if not latest_edited:
            raise ValueError(f'VC {self.id} does not have a last edited version')
        return latest_edited.publish(share_level=share_level, user=user, vc=self)

    @property
    def last_edited_version(self) -> 'ClassificationModification':
        return ClassificationModification.objects.filter(classification=self, is_last_edited=True).first()

    def latest_modification_for_user(self, user: User, exclude_withdrawn: bool = True) -> 'ClassificationModification':
        return ClassificationModification.latest_for_user(user, classification=self,
                                                          exclude_withdrawn=exclude_withdrawn).first()

    @property
    def last_published_version(self) -> 'ClassificationModification':
        return ClassificationModification.objects.filter(classification=self, is_last_published=True) \
            .select_related('classification', 'classification__lab', 'classification__lab__organization') \
            .first()

    @property
    def last_published_sync_records(self):
        records = []
        if lpv := self.last_published_version:
            records = lpv.classificationmodificationsyncrecord_set.filter(success=True)
        return records

    @staticmethod
    def validate_evidence(evidence: dict) -> List[Dict]:
        messages = []

        for key, blob in evidence.items():
            if blob is not None:  # just required to not blow up on legacy, shouldn't have Nones in here
                validations = blob.get('validation')
                if validations is not None:
                    if not isinstance(validations, list):
                        validations = [validations]  # just required to not blow up on legacy
                    for validation in validations:
                        try:
                            valid_copy = validation.copy()
                            valid_copy['key'] = key
                            messages.append(valid_copy)
                        except:
                            pass
        # make order of messages predictable
        messages.sort(key=lambda m: m.get('key') + ':' + m.get('code'))
        return messages

    def current_state_validation(self, user: User) -> ValidationMerger:
        if self.id:
            self.refresh_from_db()
        validations_received = classification_current_state_signal.send(sender=Classification, record=self, user=user)
        multi_errors = ValidationMerger.union_send_responses(validations_received)

        return multi_errors

    def has_errors(self) -> bool:
        for blob in self.evidence.values():
            if blob is not None:  # just required to not blow up on legacy
                validations = blob.get('validation')
                if validations is not None:
                    if not isinstance(validations, list):
                        validations = [validations]  # just required to not blow up on legacy
                    for validation in validations:
                        if validation.get('severity') == 'error':
                            return True
        return False

    def validate(self) -> List[Dict]:
        return Classification.validate_evidence(self.evidence)

    @staticmethod
    def flatten(record: Mapping, prefix: str = '', ignore_none_values=False, depth=0) -> dict:
        flat = {}
        for key, value in record.items():
            if ignore_none_values and value is None:
                continue
            use_key = prefix
            if use_key:
                use_key += '.' + key
            else:
                use_key = key

            if isinstance(value, Mapping):
                flat.update(Classification.flatten(value, prefix=use_key, depth=depth + 1))
            elif key == 'value':
                if isinstance(value, list):
                    if value:
                        value = ', '.join(value)
                    else:
                        value = None
                flat[prefix] = value
            else:
                flat[use_key] = value
        return flat

    def fix_permissions(self, fix_modifications=False):
        clear_permissions(self, [self.get_write_perm()])
        # labs can edit classifications
        assign_perm(self.get_write_perm(), self.lab.group, self)
        if fix_modifications:
            for cm in self.classificationmodification_set.all():
                cm.fix_permissions()

    @property
    def unique_patient_id(self):
        patient_id = self.get(SpecialEKeys.PATIENT_ID)
        if patient_id is None or self.lab is None:
            return None
        return self.lab.group_name + '/' + patient_id

    @property
    def _evidence_key_overrides_from_evidence_fields(self) -> EvidenceKeyOverrides:
        namespaces = set()
        for namespace_key in [SpecialEKeys.ASSERTION_METHOD, SpecialEKeys.ALLELE_ORIGIN]:
            for value in self.get_value_list(namespace_key):
                if namespaces_key_found := EvidenceKeyMap.cached_key(namespace_key).option_dictionary_property(
                        "namespaces").get(value):
                    namespaces |= set(namespaces_key_found)

        return EvidenceKeyOverrides(
            namespaces=namespaces
        )

    @cached_property
    def evidence_key_overrides(self) -> EvidenceKeyOverrides:
        return EvidenceKeyOverrides.merge(
            EvidenceKeyOverrides.from_dict(self.lab.organization.classification_config),
            EvidenceKeyOverrides.from_dict(self.lab.classification_config),
            self._evidence_key_overrides_from_evidence_fields
        )

    def _clear_evidence_key_cache(self):
        # ignore PyCharm's warnings https://youtrack.jetbrains.com/issue/PY-46623/Deleting-a-functools.cachedproperty-throws-a-warning-that-it-cannot-be-deleted
        clear_cached_property(self, 'evidence_key_overrides')
        clear_cached_property(self, 'evidence_keys')

    def as_json(self, params: ClassificationJsonParams) -> dict:
        from classification.models.classification_json import populate_classification_json
        return populate_classification_json(self, params)

    def lowest_share_level(self, user) -> ShareLevel:
        """ lower = less restricted """

        if user.is_superuser:
            return ShareLevel.CURRENT_USER

        is_in_lab = user.groups.filter(name=self.lab.group_name).exists()
        is_in_org = user.groups.filter(name=self.lab.organization.group_name).exists()

        if self.user == user:
            return ShareLevel.CURRENT_USER
        if is_in_lab:
            return ShareLevel.LAB
        if is_in_org:
            return ShareLevel.INSTITUTION
        if user.is_authenticated:
            return ShareLevel.ALL_USERS
        return ShareLevel.PUBLIC

    def get_visible_evidence(self, evidence, lowest_share_level: ShareLevel) -> Dict[str, Dict]:
        """ Driven by EvidenceKey.max_share_level """

        if lowest_share_level.index == 0:  # No restrictions
            return evidence

        visible_ekeys = self.evidence_keys.share_level_and_higher(lowest_share_level)
        visible_keys = {ek.key for ek in visible_ekeys}
        visible_evidence = {}
        for k, v in evidence.items():
            if k in visible_keys:
                visible_evidence[k] = v
            elif 'value' in v or 'explain' in v or 'note' in v:
                visible_evidence[k] = {'value': "(hidden)", 'hidden': True}
        return visible_evidence

    def get_allele_info_dict(self) -> Optional[Dict[str, Any]]:
        allele_info_dict = {}
        if allele_info := self.allele_info:
            resolved_dict = {
                "allele_id": allele_info.allele_id,
                "allele_info_id": allele_info.id,
                "allele_info_status": allele_info.status,
                "status": allele_info.status,
                "include": allele_info.latest_validation.include if allele_info.latest_validation else None,
                "variant_coordinate": allele_info.variant_coordinate
            }

            if (genome_build := self.get_genome_build_opt()) and \
                    (preferred_build := allele_info[genome_build]) and \
                    (c_hgvs := preferred_build.c_hgvs_obj):
                resolved_dict.update(c_hgvs.to_json())
            elif c_hgvs_raw := self.get(SpecialEKeys.C_HGVS):
                resolved_dict.update(CHGVS(c_hgvs_raw).to_json())

            include = False
            if latest_validation := allele_info.latest_validation:
                include = latest_validation.include

            resolved_dict["include"] = include
            if warning_icon := ImportedAlleleInfo.icon_for(status=allele_info.status, include=include):
                resolved_dict.update(warning_icon.as_json())

            allele_info_dict["resolved"] = resolved_dict

            genome_builds = {}
            for variant_info in allele_info.resolved_builds:
                genome_builds[variant_info.genome_build.name] = {
                    'variant_id': variant_info.variant_id,
                    SpecialEKeys.C_HGVS: variant_info.c_hgvs
                }

            if genome_builds:
                allele_info_dict["genome_builds"] = genome_builds

        return allele_info_dict

    @staticmethod
    def get_url_for_pk(pk):
        return reverse('view_classification', kwargs={'classification_id': pk})

    def get_absolute_url(self):
        return self.get_url_for_pk(self.pk)

    def get_edit_url(self) -> str:
        return self.get_absolute_url() + "?edit=true"

    @staticmethod
    def get_q_for_gene(gene: Gene) -> Q:
        match_gene = Q(**{"variant__" + VariantAnnotation.GENE_COLUMN: gene})
        match_evidence = Q(evidence__ensembl_gene_id__value=gene.pk)
        return match_gene | match_evidence

    @staticmethod
    def get_classifications_qs(user: User, clinical_significance_list: Iterable[str] = None,
                               lab_list: Iterable[Lab] = None) -> QuerySet:
        cm_qs = ClassificationModification.latest_for_user(user, published=True)
        if clinical_significance_list:
            cm_qs = cm_qs.filter(clinical_significance__in=clinical_significance_list)
        qs = Classification.objects.filter(pk__in=cm_qs.values('classification'))
        if lab_list:
            qs = qs.filter(lab__in=lab_list)
        return qs

    @staticmethod
    def get_variant_q(user: User, genome_build: GenomeBuild,
                      clinical_significance_list: Iterable[str] = None,
                      lab_list: Iterable[Lab] = None) -> Q:
        """ returns a Q object filtering variants to those with a PUBLISHED classification
            (optionally classification in clinical_significance_list """
        vc_qs = Classification.get_classifications_qs(user, clinical_significance_list, lab_list)
        return Classification.get_variant_q_from_classification_qs(vc_qs, genome_build)

    @staticmethod
    def get_variant_q_from_classification_qs(vc_qs, genome_build: GenomeBuild) -> Q:
        va_qs = VariantAllele.objects.filter(genome_build=genome_build,
                                             allele__in=vc_qs.values_list("allele"))
        variant_ids = va_qs.values_list("variant_id", flat=True)
        return Q(id__in=list(variant_ids))  # List is much faster than inner query...

    @staticmethod
    def annotate_with_variant_sort(classifications_qs: QuerySet, genome_build: GenomeBuild, name="variant_sort"):
        """ Annotate Classification queryset you can use to order by genome build position """

        variant_qs: QuerySet
        if classifications_qs.model == Classification:
            variant_qs = Variant.objects.filter(variantallele__allele=OuterRef("variant__variantallele__allele"),
                                                variantallele__genome_build=genome_build)
        else:
            variant_qs = Variant.objects.filter(
                variantallele__allele=OuterRef("classification__variant__variantallele__allele"),
                variantallele__genome_build=genome_build)

        variant_qs = variant_qs.annotate(
            padded_contig=LPad("locus__contig__name", 2, Value("0")),
            padded_position=LPad(Cast("locus__position", output_field=TextField()), 9, Value("0")),
            variant_sort=Concat("padded_contig", "padded_position"),
        )

        return classifications_qs.annotate(**{name: Subquery(variant_qs[:1].values_list("variant_sort"))})

    def get_other_classifications_summary_for_variant(self, user: User) -> str:
        other_classifications_summary = None
        if self.variant and self.variant.allele:
            latest_others_for_variant = ClassificationModification.latest_for_user(user=user,
                                                                                   allele=self.variant.allele,
                                                                                   published=True)
            latest_others_for_variant = latest_others_for_variant.exclude(classification__id=self.id)

            counts = Counter(
                [vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE, "Unclassified") for vcm in latest_others_for_variant])
            classifications = []
            for counted in counts.most_common():
                count = counted[1]
                if count:
                    cs = counted[0]
                    cs_label = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(cs)
                    classifications.append(f"{cs_label} x{count}")
            other_classifications_summary = ", ".join(classifications)

        return other_classifications_summary

    def get_variant_for_build(self, genome_build: GenomeBuild) -> Optional[Variant]:
        if self.variant and self.variant.allele:
            return self.variant.allele.variant_for_build_optional(genome_build)
        return None

    def get_variant_annotation(self, variant_annotation_version: VariantAnnotationVersion) -> Optional[
        VariantAnnotation]:
        if variant := self.get_variant_for_build(variant_annotation_version.genome_build):
            return variant.variantannotation_set.filter(version=variant_annotation_version).first()

    def c_hgvs_all(self) -> List[CHGVS]:
        all_chgvs: List[CHGVS] = []
        for genome_build in GenomeBuild.builds_with_annotation_cached():
            if text := self.get_c_hgvs(genome_build):
                chgvs = CHGVS(full_c_hgvs=text)
                chgvs.genome_build = genome_build
                chgvs.is_normalised = True
                all_chgvs.append(chgvs)
        return all_chgvs

    def c_hgvs_best(self, preferred_genome_build: GenomeBuild) -> CHGVS:
        if c_hgvs_str := self.get_c_hgvs(preferred_genome_build):
            c_hgvs = CHGVS(c_hgvs_str)
            c_hgvs.genome_build = preferred_genome_build
            c_hgvs.is_normalised = True
            c_hgvs.is_desired_build = True
            return c_hgvs
        for alt_genome_build in GenomeBuild.builds_with_annotation_cached():
            if preferred_genome_build == alt_genome_build:
                continue
            if c_hgvs_str := self.get_c_hgvs(alt_genome_build):
                c_hgvs = CHGVS(c_hgvs_str)
                c_hgvs.genome_build = alt_genome_build
                c_hgvs.is_normalised = True
                c_hgvs.is_desired_build = False
                return c_hgvs
        c_hgvs = CHGVS(self.get(SpecialEKeys.C_HGVS) or "")
        try:
            c_hgvs.genome_build = self.get_genome_build()
        except ValueError:
            pass
        c_hgvs.is_normalised = False
        c_hgvs.is_desired_build = preferred_genome_build == c_hgvs.genome_build
        return c_hgvs

    def _generate_c_hgvs(self, genome_build: GenomeBuild) -> str:
        variant = self.get_variant_for_build(genome_build)
        hgvs_matcher = HGVSMatcher(genome_build=genome_build)

        c_hgvs: str = None
        if variant:
            transcript_id = None
            try:
                transcript_id = self.transcript
                hgvs_variant = hgvs_matcher.variant_to_hgvs_variant(variant, transcript_id)
                c_hgvs = hgvs_variant.format()
            except Exception:
                # can't map between builds
                report_exc_info(extra_data={
                    "genome_build": genome_build.name,
                    "variant": str(variant),
                    "transcript_id": transcript_id
                })
        return c_hgvs

    def get_c_hgvs(self, genome_build: GenomeBuild, use_compat: bool = False) -> Optional[str]:
        if genome_build == genome_build.grch37():
            return self.chgvs_grch37 if not use_compat else self.chgvs_grch37_compat
        if genome_build == genome_build.grch38():
            return self.chgvs_grch38 if not use_compat else self.chgvs_grch38_compat
        return self._generate_c_hgvs(genome_build)

    def __str__(self) -> str:
        parts = [f"({str(self.id)})"]
        genome_build = GenomeBuildManager.get_current_genome_build()
        cached_c_hgvs = self.get_c_hgvs(genome_build=genome_build)
        if not cached_c_hgvs:
            cached_c_hgvs = self.get(SpecialEKeys.C_HGVS)
        parts.append(cached_c_hgvs or "No c.HGVS")

        clinical_significance = self.get_clinical_significance_display() or "No Data"
        parts.append(clinical_significance)

        if self.withdrawn:
            parts.append("WITHDRAWN")

        return " ".join(parts)

    @staticmethod
    def check_can_create_no_classification_via_web_form(_user: User):
        if not settings.CLASSIFICATION_WEB_FORM_CREATE_ALLOW_NO_VARIANT:
            raise CreateNoClassificationForbidden()


class ClassificationModification(GuardianPermissionsMixin, EvidenceMixin, models.Model):
    classification = models.ForeignKey(Classification, on_delete=CASCADE)
    user = models.ForeignKey(User, on_delete=PROTECT)  # One who did last change, may not be classification.user
    created = DateTimeUTCField(db_index=True, auto_now_add=True)
    modified = DateTimeUTCField(auto_now=True)
    source = models.TextField(null=False)

    @property
    def source_enum(self) -> SubmissionSource:
        return SubmissionSource(self.source) if self.source else SubmissionSource.FORM

    delta = models.JSONField(null=False, blank=True, default=dict)
    published = models.BooleanField(null=False, default=False)
    published_evidence = models.JSONField(null=True, blank=True, default=None)
    # changing the share level does not change the actual logic of sharing
    # it is useful if we have to rework permissions though
    share_level = models.CharField(max_length=16, choices=ShareLevel.choices(), null=True, blank=True)

    @property
    def imported_c_hgvs_obj(self) -> CHGVS:
        if c_hgvs := self.get(SpecialEKeys.C_HGVS):
            # remove any white space inside the c.HGVS
            c_hgvs = re.sub(r'\s+', '', c_hgvs)
            c_hgvs_obj = CHGVS(c_hgvs)
            try:
                c_hgvs_obj.genome_build = self.get_genome_build()
            except ValueError:
                pass
            return c_hgvs_obj

    @property
    def allele_origin_bucket_obj(self) -> AlleleOriginBucket:
        return self.classification.allele_origin_bucket

    @property
    def condition_text(self):
        if crd := self.classification.condition_resolution_dict:
            return crd.get('display_text')
        return self.get(SpecialEKeys.CONDITION)

    @property
    def condition_resolution_dict_fallback(self) -> ConditionResolvedDict:
        if resolved := self.classification.condition_resolution:
            return resolved
        return ConditionResolvedDict(
            display_text=self.get(SpecialEKeys.CONDITION)
        )

    @cached_property
    def condition_resolution_obj_fallback(self) -> ConditionResolved:
        # Would rather replace condition_resolution_obj with this
        # but don't know if that would break anything
        if cr_obj := self.classification.condition_resolution_obj:
            return cr_obj
        else:
            return ConditionResolved(terms=[], plain_text=self.get(SpecialEKeys.CONDITION))

    @staticmethod
    def column_name_for_build(genome_build: GenomeBuild, suffix: str = 'c_hgvs'):
        return ImportedAlleleInfo.column_name_for_build(genome_build, "classification__allele_info", suffix)

    @property
    def share_level_enum(self) -> ShareLevel:
        if self.share_level:
            return ShareLevel(self.share_level)

        # for historical data before we tracked share_level
        for sl in [ShareLevel.PUBLIC,
                   ShareLevel.ALL_USERS,
                   ShareLevel.INSTITUTION,
                   ShareLevel.LAB]:
            if self.can_view(sl.group(lab=self.classification.lab)):
                return sl
        return ShareLevel.CURRENT_USER

    @share_level_enum.setter
    def share_level_enum(self, value: ShareLevel):
        self.share_level = value.value

    clinical_significance = models.CharField(max_length=1, choices=ClinicalSignificance.CHOICES, null=True, blank=True)

    is_last_published = models.BooleanField(db_index=True, null=False, blank=True, default=False)
    is_last_edited = models.BooleanField(db_index=True, null=False, blank=True, default=False)

    @classmethod
    def order_by_evidence(cls, key_id: str) -> RawSQL:
        return RawSQL('cast(published_evidence->>%s as jsonb)->>%s', (key_id, 'value'))

    def __str__(self) -> str:
        return self.id_str

    @property
    def _evidence(self):
        return self.evidence

    @property
    def id_str(self) -> str:
        return self.classification.id_str + '.' + str(self.created.timestamp())

    @property
    def cr_lab_id(self) -> str:
        if settings.CLASSIFICATION_ID_OVERRIDE_PREFIX:
            return f"CR_{self.classification.id_str}"
        return self.classification.lab_record_id

    @property
    def curated_date(self) -> datetime:
        return CuratedDate(self).date

    @cached_property
    def curated_date_check(self) -> 'CuratedDate':
        return CuratedDate(self)

    @property
    def lab(self) -> Lab:
        return self.classification.lab

    def as_json(self,
                params: ClassificationJsonParams) -> dict:

        params.version = self
        return self.classification.as_json(params)

    def get_absolute_url(self):
        return reverse('view_classification',
                       kwargs={'classification_id': str(self.classification.id) + '.' + str(self.created.timestamp())})

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cached_evidence = None

    @classmethod
    def latest_for_user(cls,
                        user: User,
                        classification: Optional[Classification] = None,
                        variant: Optional[Union[Variant, Iterable[Variant]]] = None,
                        allele: Allele = None,
                        published: bool = None,
                        clinical_context: Optional['ClinicalContext'] = None,
                        exclude_withdrawn: bool = True,
                        shared_only: bool = False,
                        allele_origin_bucket: Optional[AlleleOriginBucket] = None,
                        allele_origin_buckets: Optional[set[AlleleOriginBucket]] = None,
                        exclude_external_labs: bool = False,
                        **kwargs) -> QuerySet['ClassificationModification']:
        """

        :param user: The active user
        :param classification: If provided will only get modifications for this classification
        :param variant: WARNING - use allele, If provided, will only get classifications linked to this variant
        :param allele: all modifications will be ones linked to this allele
        :param published: If true, only the last published record for a classification will be returned
        :param clinical_context: If provided, only modifications linked to the clinical context will be returned
        :param exclude_withdrawn: True by default, if false will include withdrawn records
        :param shared_only: False by default, if true will only include records at discordant levels
        :param allele_origin_bucket: None by default, if Germline/Somatic/Unknown will only return classification from those buckets
        :param allele_origin_buckets: If we want multiple buckets, provide a set of them
        :param exclude_external_labs: If True, will not consider external labs to select from
        :param kwargs: Any other parameters will be used for filtering
        :return:
        """
        qs: QuerySet[ClassificationModification] = get_objects_for_user(user, cls.get_read_perm(), klass=cls,
                                                                        accept_global_perms=True)
        qs = qs.order_by('-created', 'classification')
        if kwargs:
            qs = qs.filter(**kwargs)

        if classification:
            qs = qs.filter(classification=classification)
        if variant is not None:
            if isinstance(variant, (Variant, int)):
                qs = qs.filter(classification__variant=variant)
            else:
                qs = qs.filter(classification__variant__in=variant)
        if allele:
            qs = qs.filter(classification__variant__in=allele.variants)
        if clinical_context:
            qs = qs.filter(classification__clinical_context=clinical_context)
        if exclude_withdrawn:
            qs = qs.filter(classification__withdrawn=False)
        if shared_only:
            qs = qs.filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS)
        if allele_origin_bucket:
            qs = qs.filter(classification__allele_origin_bucket=allele_origin_bucket)
        elif allele_origin_buckets:
            qs = qs.filter(classification__allele_origin_bucket__in=allele_origin_buckets)
        if exclude_external_labs:
            qs = qs.filter(classification__lab__external=False)

        if published:
            qs = qs.select_related(
                'classification',
                'classification__lab__organization',
                'classification__allele_info__grch37',
                'classification__allele_info__grch38',
                'classification__user'
            )
            qs = qs.filter(is_last_published=True)
        else:
            if classification:  # if we're only looking at one classification, grab the latest the user can see
                # not guaranteed to be is_last_edited
                if latest_single := qs.order_by('-created').first():
                    qs = ClassificationModification.objects.filter(pk=latest_single.pk)
                else:
                    return ClassificationModification.objects.none()

            else:
                # FIXME this will not get classifications with outstanding edits that the user can see
                qs = qs.filter(is_last_edited=True)

        return qs

    def is_edit_appendable(self, user: User, source: SubmissionSource) -> bool:
        """
        Return if we can edit this record instead of creating a new one
        """
        edit_window = datetime.utcnow().replace(tzinfo=timezone.utc) - timedelta(minutes=1)

        return \
                self.source_enum == source \
                and self.created >= edit_window \
                and self.user == user \
                and not self.published

    @transaction.atomic()
    def publish(self, share_level: Union[ShareLevel, str, int], user: User, vc: 'Classification') -> bool:
        """
        :param share_level: The share level we want to publish as
        :param user: The user who initiated the publishing
        :param vc: The variant classification we're publishing - even though it's redundant, but we get the same object
        instance and don't reload it from the database
        :return a boolean indicating if a new version was published (false if no change was required)
        """
        old_share_level = vc.share_level_enum
        share_level = ShareLevel.from_key(share_level)
        if old_share_level == share_level and self.is_last_published:
            # nothing to do, do not re-publish
            return False

        lab = self.classification.lab
        vc.share_level_enum = share_level

        # get the new share level (have to pass in lab to determine institution
        # if share level happens to be that
        group = share_level.group(lab=lab)

        previously_published = vc.last_published_version

        self.published = True
        self.published_evidence = vc.evidence
        self.share_level_enum = share_level

        if previously_published and previously_published.id != self.id:
            # make sure to unmark previous record as the last published
            previously_published.is_last_published = False
            previously_published.save(update_fields=['is_last_published'])

        self.is_last_published = True
        self.save()

        if group:
            assign_perm(self.get_read_perm(), group, self)

        vc.summary = ClassificationSummaryCalculator(self).cache_dict()
        vc.save()

        get_timer().tick("Published Modification")

        classification_post_publish_signal.send(
            sender=Classification,
            classification=vc,
            previously_published=previously_published,
            previous_share_level=old_share_level,
            newly_published=self,
            user=user)

        vc.refresh_from_db()
        return True

    def fix_permissions(self):
        clear_permissions(self, [self.get_read_perm()])
        assign_perm(self.get_read_perm(), self.classification.lab.group, self)
        if self.published:
            if group := self.share_level_enum.group(lab=self.classification.lab):
                assign_perm(self.get_read_perm(), group, self)

    @cached_property
    def previous(self) -> Optional['ClassificationModification']:
        return ClassificationModification.objects.filter(classification=self.classification,
                                                         created__lt=self.created).order_by('-created').first()

    @property
    def index(self):
        return ClassificationModification.objects.filter(classification=self.classification,
                                                         created__lt=self.created).count()

    @property
    def evidence(self) -> Dict[str, Dict]:
        if self.cached_evidence is None:
            if self.published_evidence is not None:
                self.cached_evidence = self.published_evidence
            else:
                data = {}
                patches_since_evidence = []

                for vcm in ClassificationModification.objects.filter(classification=self.classification,
                                                                     created__lte=self.created) \
                        .order_by('-created').values('delta', 'published_evidence'):
                    published_evidence = vcm['published_evidence']
                    if published_evidence is not None:
                        data = published_evidence
                        break
                    patches_since_evidence.append(vcm['delta'])

                patches_since_evidence.reverse()
                for patch in patches_since_evidence:
                    Classification.patch_with(target=data, patch=patch, tidy_nones=True)
                self.cached_evidence = data
        return self.cached_evidence

    def get_visible_evidence(self, user: User) -> Dict[str, Dict]:
        """ Driven by EvidenceKey.max_share_level """
        lowest_share_level = self.classification.lowest_share_level(user)
        return self.classification.get_visible_evidence(self.evidence, lowest_share_level)

    def c_hgvs_best(self, genome_build: GenomeBuild) -> CHGVS:
        return self.classification.c_hgvs_best(preferred_genome_build=genome_build)

    def is_significantly_equal(self, other: 'ClassificationModification', care_about_explains: bool = False) -> bool:
        """
        Determines as far as a user is concerned if two versions of a classification are equal based on evidence key
        data and share level.
        Would typically be used to determine if we should show the column in history.
        Raises an error if run on different classification ids.
        :param other: Another modification
        :param care_about_explains: True = change in only 'explain' should cause records not to be considered the same
        :return: True if the user should see these as the same, False otherwise
        """
        if self.classification_id != other.classification_id:
            raise ValueError("Called is_significantly_equal on modifications for different classifications")

        if self.share_level_enum != other.share_level_enum:
            return False

        def is_prop_equal(a: Any, b: Any) -> bool:
            if a == b:
                return True
            if not isinstance(a, list) and not isinstance(b, list):
                return False

            # one or both of these values is an array, check array equality
            # useful for if a key changed from a select to a multi-select
            if not isinstance(a, list):
                if isinstance(a, str):
                    a = [s.strip() for s in a.split(',')]
                else:
                    a = [a]
            if not isinstance(b, list):
                if isinstance(b, str):
                    b = [s.strip() for s in b.split(',')]
                else:
                    b = [b]
            if len(a) != len(b):
                return False
            a.sort()
            b.sort()
            for a_sub, b_sub in zip(a, b):
                if a_sub != b_sub:
                    return False
            return True

        # we only care about changes to actual keys
        e_keys = EvidenceKeyMap.instance()
        self_evidence = self.evidence
        other_evidence = other.evidence

        check_keys = set(self_evidence.keys()).union(other_evidence.keys())
        # the following keys do not represent data and should not be considered
        # changes worth displaying
        check_keys.discard(SpecialEKeys.OWNER)
        check_keys.discard(SpecialEKeys.SOURCE_ID)

        check_props = ['value', 'note']
        if care_about_explains:
            check_props.append('explain')

        for key in check_keys:
            # ignore invalid key values
            if e_keys.get(key).is_dummy:
                continue
            source_blob = self_evidence.get(key) or {}
            other_blob = other_evidence.get(key) or {}
            if source_blob != other_blob:
                for attribute in check_props:
                    if not is_prop_equal(a=source_blob.get(attribute), b=other_blob.get(attribute)):
                        return False
        return True


@dataclass
class ClassificationConsensus:
    modification: ClassificationModification
    label: str = "no-label"
    default_suggestion: bool = False

    @staticmethod
    def all_consensus_candidates(allele: Allele, user: User) -> ['ClassificationConsensus']:
        us = UserSettingsManager.get_user_settings(user)

        default_allele_origin_filter = us.allele_origin_focus
        default_allele_origin_bucket: Optional[AlleleOriginBucket] = None
        if default_allele_origin_filter.value in (AlleleOriginBucket.GERMLINE.value, AlleleOriginBucket.SOMATIC.value):
            default_allele_origin_bucket = AlleleOriginBucket(default_allele_origin_filter.value)

        results = []
        for identifier, allele_origin_bucket in [("Latest Germline", AlleleOriginBucket.GERMLINE),
                                                 ("Latest Somatic", AlleleOriginBucket.SOMATIC)]:
            if candidates := list(ClassificationModification.latest_for_user(user=user, allele=allele, published=True,
                                                                             exclude_external_labs=True,
                                                                             allele_origin_bucket=allele_origin_bucket).all()):
                candidates.sort(key=lambda vcm: vcm.curated_date_check, reverse=True)
                default_suggestion = allele_origin_bucket == default_allele_origin_bucket
                results.append(ClassificationConsensus(modification=candidates[0], label=identifier,
                                                       default_suggestion=default_suggestion))
        return results

    @cached_property
    def consensus_patch(self) -> VCPatch:
        keys = EvidenceKeyMap.cached()

        evidence = self.modification.published_evidence
        consensus: Dict[str, Any] = {}

        # default allele origin - don't use copy consensus because that would copy "likely somatic" etc
        if allele_origin_bucket := self.modification.classification.allele_origin_bucket:
            allele_origin = None
            if allele_origin_bucket == AlleleOriginBucket.GERMLINE:
                allele_origin = "germline"
            elif allele_origin_bucket == AlleleOriginBucket.SOMATIC:
                allele_origin = "somatic"

            if allele_origin:
                consensus["allele_origin.value"] = allele_origin

        for key in (ekey.key for ekey in keys.all_keys if ekey.copy_consensus):
            for part in ['value', 'note']:
                blob = evidence.get(key)
                if blob is None:
                    continue
                part_value = blob.get(part)
                if part_value is None or (isinstance(part_value, list) and not part_value):
                    continue
                consensus[key + '.' + part] = part_value

        patch: Dict[str, Dict[str, Any]] = nest_dict(consensus)
        return patch


class ClassificationDateType(StrEnum):
    CURATION = ""  # default
    VERIFIED = "Verified"
    SAMPLE_DATE = "Sample Date"
    CREATED = "Created"


@dataclass(frozen=True)
class ClassificationDate:
    date_type: ClassificationDateType
    datetime: Optional[datetime] = None
    date: Optional[date] = None

    def __post_init__(self):
        if not self.date and not self.datetime:
            raise ValueError("Either datetime or date must be provided")

    @property
    def name(self):
        return self.date_type.value

    @property
    def value(self) -> Union[date, datetime]:
        if self.date:
            return self.date
        return self.datetime

    @property
    def date_str(self) -> str:
        return Classification.to_date_str(self.value)

    def __lt__(self, other: 'ClassificationDate'):
        if self.datetime and other.datetime:
            return self.datetime < other.datetime
        else:
            return (self.date or date.min) < (other.date or date.min)


CLASSIFICATION_DATE_REGEX = re.compile(r"(?P<year>[0-9]{4})-(?P<month>[0-9]{2})-(?P<day>[0-9]{2})")


class CuratedDate:
    """
    CuratedDate is a bit misleading, keeps the most relevant filled in date for a classification.
    Also allows two classification modifications can be compared, and always having a winner
    (subsequent less important dates are fallen back to if other dates are equal)
    """

    def __init__(self, modification: ClassificationModification):
        self._modification = modification

    @cached_property
    def timezone(self):
        return gettz(settings.TIME_ZONE)

    def convert_date(self, evidence_key) -> Optional[date]:
        if date_str := self._modification.get(evidence_key):
            if m := CLASSIFICATION_DATE_REGEX.match(date_str):
                try:
                    return date(year=int(m.group("year")), month=int(m.group("month")), day=int(m.group("day")))
                except ValueError:
                    # an invalid date should already cause a warning on the classification form
                    pass

    @cached_property
    def curation_date(self) -> Optional[ClassificationDate]:
        if date_val := self.convert_date(SpecialEKeys.CURATION_DATE):
            return ClassificationDate(date_type=ClassificationDateType.CURATION, date=date_val)

    @cached_property
    def curated_verified_date(self) -> Optional[ClassificationDate]:
        if date_val := self.convert_date(SpecialEKeys.CURATION_VERIFIED_DATE):
            return ClassificationDate(date_type=ClassificationDateType.VERIFIED, date=date_val)

    @cached_property
    def sample_date(self) -> Optional[ClassificationDate]:
        if date_val := self.convert_date(SpecialEKeys.SAMPLE_DATE):
            return ClassificationDate(date_type=ClassificationDateType.SAMPLE_DATE, date=date_val)

    @cached_property
    def created_date(self) -> ClassificationDate:
        return ClassificationDate(date_type=ClassificationDateType.CREATED, datetime=self._modification.classification.created.astimezone(self.timezone))

    @cached_property
    def relevant_date(self) -> ClassificationDate:
        """
        Return the most relevant date (for reporting) along with a name for the date if it isn't the Sample Date
        """
        if curation_date := self.curation_date:
            return curation_date
        if curation_verified_date := self.curated_verified_date:
            return curation_verified_date
        if sample_date := self.sample_date:
            return sample_date
        return self.created_date

    @property
    def name(self) -> str:
        return self.relevant_date.name

    @property
    def value(self) -> Union[date, datetime]:
        return self.relevant_date.value

    @property
    def date(self) -> Union[date, datetime]:
        return self.value

    def __lt__(self, other: 'CuratedDate') -> bool:
        my_relevant_date = self.relevant_date
        other_relevant_date = other.relevant_date

        if my_relevant_date.value != other_relevant_date.value:
            return my_relevant_date < other_relevant_date

        # Most relevant dates are the same, now go down the chain until we find a difference
        def direction(date_1: ClassificationDate, date_2: ClassificationDate):
            if date_1 is None and date_2 is not None:
                return -1
            elif date_1 is not None and date_2 is None:
                return 1
            if date_1 is None and date_2 is None:
                return 0

            if date_1.value == date_2.value:
                return 0
            elif date_1 < date_2:
                return -1
            else:
                return 1

        if diff := direction(self.curation_date, other.curation_date):
            return diff < 0
        if diff := direction(self.curated_verified_date, other.curated_verified_date):
            return diff < 0
        if diff := direction(self.sample_date, other.sample_date):
            return diff < 0
        if diff := direction(self.created_date, other.created_date):
            return diff < 0

        # should not happen as two records should never have the same created date
        # but just in case
        return self._modification.classification_id < other._modification.classification_id
