import copy
import re
import uuid
from collections import Counter, namedtuple
from dataclasses import dataclass, field
from datetime import datetime, timezone, timedelta
from enum import Enum
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
    CRITERIA_NOT_MET, ValidationCode, CriteriaEvaluation
from classification.models.classification_import_run import ClassificationImportRun
from classification.models.classification_patcher import patch_fuzzy_age
from classification.models.classification_utils import \
    ValidationMerger, ClassificationJsonParams, PatchMeta, ClassificationPatchResponse
from classification.models.classification_variant_info_models import ImportedAlleleInfo, ImportedAlleleInfoStatus
from classification.models.evidence_key import EvidenceKeyValueType, \
    EvidenceKey, EvidenceKeyMap, VCDataDict, WipeMode, VCDataCell
from classification.models.evidence_mixin import EvidenceMixin, VCPatch
from classification.models.flag_types import classification_flag_types
from flags.models import Flag, FlagPermissionLevel, FlagStatus
from flags.models.models import FlagsMixin, FlagCollection, FlagTypeContext, \
    flag_collection_extra_info_signal, FlagInfos
from genes.hgvs import HGVSMatcher, CHGVS, HGVSNameExtra
from genes.models import Gene
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.guardian_utils import clear_permissions
from library.log_utils import report_exc_info, report_event
from library.utils import empty_to_none, nest_dict, cautious_attempt_html_to_text, DebugTimer, \
    invalidate_cached_property, md5sum_str
from ontology.models import OntologyTerm, OntologySnake, OntologyTermRelation
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Variant, Lab, Sample
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import AlleleSource, Allele, VariantAllele

ChgvsKey = namedtuple('CHGVS', ['short', 'column', 'build'])

classification_validation_signal = django.dispatch.Signal()  # args: "classification", "patch_meta", "key_map"
classification_current_state_signal = django.dispatch.Signal()  # args: "user"
classification_post_publish_signal = django.dispatch.Signal()  # args: "classification", "previously_published", "previous_share_level", "newly_published", "user", "debug_timer"
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
    e.g. referring to a non existent lab
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
        return Variant.objects.filter(pk__in=ImportedAlleleInfo.objects.filter(classification_import=self).values_list('matched_variant', flat=True))

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
    """ Used to reload all Classifications (for upgrades etc) """

    script = models.TextField()  # Set if run from script, can check hasn't been re-run
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    git_hash = models.TextField()

    def get_genome_build(self) -> GenomeBuild:
        return self.genome_build

    def get_variants_qs(self) -> QuerySet[Variant]:
        # Note: This deliberately only gets classifications where the submitting variant was against this genome build
        # ie we don't use Classification.get_variant_q_from_classification_qs() to get liftovers
        contigs_q = Variant.get_contigs_q(self.genome_build)

        not_lifted_over_variant_ids = ImportedAlleleInfo.objects.filter(matched_variant__isnull=False, allele__isnull=True).values_list('allele_id', flat=True)

        return Variant.objects.filter(contigs_q, id__in=not_lifted_over_variant_ids)

    def liftover_complete(self, genome_build: GenomeBuild):
        ImportedAlleleInfo.relink_variants()


@receiver(flag_collection_extra_info_signal, sender=FlagCollection)
def get_extra_info(flag_infos: FlagInfos, user: User, **kwargs) -> None:  # pylint: disable=unused-argument
    """
    Allows us to provide extra info for FlagCollections attached to Classification
    e.g. linking to the appropriate allele page, discordance report etc
    :param flag_infos: Information on the flag collections being displayed to the user.
    Populates this with the extra info
    :param user: The current user
    :param kwargs: Required by @receiver
    """
    from classification.models.discordance_models import DiscordanceReportClassification
    from classification.enums.discordance_enums import DiscordanceReportResolution

    vcs = Classification.objects.filter(flag_collection__in=flag_infos.ids).select_related('lab')
    drcs = DiscordanceReportClassification.objects.filter(classification_original__classification__in=vcs,
                                                          report__resolution=DiscordanceReportResolution.ONGOING) \
        .values_list('classification_original__classification', 'report')
    drcs_dict = {}
    for drc in drcs:
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
    resolved_terms: List[ConditionResolvedTermDict]
    resolved_join: str


@dataclass(frozen=True)
class ConditionResolved:
    terms: List[OntologyTerm]
    join: Optional['MultiCondition'] = None
    plain_text: Optional[str] = field(default=None)  # fallback, not populated in all contexts

    @property
    def summary(self) -> str:
        text = ", ".join([term.id for term in self.terms])
        if join := self.join:
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
        return ConditionResolved(terms=terms, join=join)

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
                    descendant_relationships = OntologySnake.check_if_ancestor(descendant=self_mondo, ancestor=other_mondo)
                    return bool(descendant_relationships)

            # terms cant be converted to MONDO and not exact match, just return False
            return False

    @staticmethod
    def more_general_term_if_related(resolved_1: 'ConditionResolved', resolved_2: 'ConditionResolved') -> Optional['ConditionResolved']:
        more_general: Optional[ConditionResolved] = None
        if resolved_1.is_same_or_more_specific(resolved_2):
            more_general = resolved_2
        elif resolved_2.is_same_or_more_specific(resolved_1):
            more_general = resolved_1

        if more_general:
            # if presented with different types, and we can switch over to MONDO, do so
            if not more_general.is_multi_condition and resolved_1.single_term.ontology_service != resolved_2.single_term.ontology_service:
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

    def to_json(self) -> ConditionResolvedDict:
        jsoned: ConditionResolvedDict
        if self.terms:
            from classification.models import MultiCondition

            def format_term(term: OntologyTerm) -> str:
                if name := term.name:
                    return f"{term.id} {name}"
                return term.id

            terms = self.terms
            text = ", ".join([format_term(term) for term in terms])
            sort_text = ", ".join([term.name for term in terms]).lower()
            join: Optional[MultiCondition] = None
            if len(terms) > 1:
                join = self.join or MultiCondition.NOT_DECIDED
                text = f"{text}; {self.join.label}"

            resolved_term_dicts: List[ConditionResolvedTermDict] = [ConditionResolved.term_to_dict(term) for term in self.terms]

            jsoned: ConditionResolvedDict = {
                "resolved_terms": resolved_term_dicts,
                "resolved_join": join,
                "display_text": text,
                "sort_text": sort_text
            }
            return jsoned
        else:
            jsoned: ConditionResolvedDict = {
                "display_text": self.plain_text,
                "sort_text": self.plain_text
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


class Classification(GuardianPermissionsMixin, FlagsMixin, EvidenceMixin, TimeStampedModel):
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
    """ Keeps links to common builds (37, 38) for quick access to c.hgvs, transcript etc. Is shared between classifications with same import data """

    @property
    def allele_object(self) -> Allele:
        """ The new preferred way to reference the allele, so we can eventually remove allele from the classification object """
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
    """ If created from a variant and auto-populated, with which version of annotations. If null was created via import """

    clinical_context = models.ForeignKey('ClinicalContext', null=True, blank=True, on_delete=SET_NULL)
    """ After being matched to a variant, this will be set to the default clinical_context for the allele
    but it can be changed to be another clinical_context (for the same allele) """

    lab_record_id = models.TextField(blank=True, null=True)
    """ Should be unique together with lab """

    evidence = models.JSONField(null=False, blank=True, default=dict)
    """ The latest evidence (should always match the content of the latest ClassificationModification.evidence) """

    withdrawn = models.BooleanField(default=False)
    """ Soft delete, if withdrawn classification should not appear in most places """

    clinical_significance = models.CharField(max_length=1, choices=ClinicalSignificance.CHOICES, null=True, blank=True)
    """ Used as an optimisation for queries """

    condition_resolution = models.JSONField(null=True, blank=True)  # of type ConditionProcessedDict

    last_source_id = models.TextField(blank=True, null=True)
    last_import_run = models.ForeignKey(ClassificationImportRun, null=True, blank=True, on_delete=SET_NULL)

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
    def chgvs_grch37_full(self) -> Optional[str]:
        try:
            return self.allele_info.grch37.c_hgvs_full
        except AttributeError:
            return None

    @property
    def chgvs_grch38_full(self) -> Optional[str]:
        try:
            return self.allele_info.grch38.c_hgvs_full
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
            return cc.name
        return 'not-matched'

    @property
    def condition_text_record(self) -> 'ConditionText':
        if ctm := self.conditiontextmatch:
            return ctm.condition_text
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

    class Meta:
        unique_together = ('lab', 'lab_record_id')

    @staticmethod
    def can_create_via_web_form(user: User):
        can_create = settings.VARIANT_CLASSIFICATION_WEB_FORM_CREATE
        return can_create and (user.is_superuser or settings.VARIANT_CLASSIFICATION_WEB_FORM_CREATE_BY_NON_ADMIN)

    @staticmethod
    def dashboard_report_new_classifications(since) -> int:
        return Classification.objects.filter(created__gte=since).count()

    @staticmethod
    def dashboard_total_shared_classifications() -> int:
        return Classification.objects.filter(lab__external=False, share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS, withdrawn=False).exclude(lab__name__icontains='legacy').count()

    @staticmethod
    def dashboard_total_unshared_classifications() -> int:
        return Classification.objects.filter(lab__external=False, withdrawn=False).exclude(lab__name__icontains='legacy').exclude(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).count()

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
        Used for the admin screen so we can show variant coordinate in listing
        """
        return self.get(SpecialEKeys.VARIANT_COORDINATE)

    @property
    def imported_genome_build(self):
        return self.get(SpecialEKeys.GENOME_BUILD)

    @property
    def imported_c_hgvs(self):
        return self.get(SpecialEKeys.C_HGVS)

    @property
    def imported_g_hgvs(self):
        return self.get(SpecialEKeys.G_HGVS)

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
    def friendly_label(self):
        return self.lab.name + ' / ' + self.lab_record_id

    @staticmethod
    def to_date(date_str: str) -> Optional[datetime]:
        if date_str:
            return datetime.strptime(date_str, '%Y-%m-%d')
        return None

    def set_withdrawn(self, user: User, withdraw: bool = False) -> bool:
        if not self.id and withdraw:
            raise ValueError('Cannot withdrawn new classification record - use delete instead')

        if self.withdrawn == withdraw:
            return False  # no change

        self.withdrawn = withdraw
        self.save()
        if withdraw:
            self.flag_collection_safe.get_or_create_open_flag_of_type(
                flag_type=classification_flag_types.classification_withdrawn,
                user=user,
                permission_check=False,
                reopen=True
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

    def ensure_allele_info_with_created(self, force_allele_info_update_check: bool = False) -> Tuple[Optional[ImportedAlleleInfo], bool]:
        created = False
        if not self.allele_info or force_allele_info_update_check:
            try:
                genome_build_patch_version = self.get_genome_build_patch_version()
            except Exception:
                # no allele info if we can't derive a genome build
                return None, False

            fields = {"imported_genome_build_patch_version": genome_build_patch_version}
            if c_hgvs := self.imported_c_hgvs:
                fields["imported_c_hgvs"] = c_hgvs
                fields["imported_md5_hash"] = md5sum_str(c_hgvs)
            elif g_hgvs := self.imported_g_hgvs:
                fields["imported_g_hgvs"] = g_hgvs
                fields["imported_md5_hash"] = md5sum_str(g_hgvs)
                fields["imported_transcript"] = self.transcript
            else:
                # need either c.hgvs, g.hgvs (in future, re-support variant_coordinate)
                return None, False

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
            allele_info.update_variant_coordinate()  # only need this for systems that were migrated when half of the AlleleInfo was done
            allele_info.update_status()
            allele_info.apply_validation()
            allele_info.save()
            if not force_update and self.allele_info.status == ImportedAlleleInfoStatus.MATCHED_ALL_BUILDS:
                return
            if matched_variant := self.variant:
                allele_info.set_variant_and_save(matched_variant, force_update=force_update)
            return True
        else:
            return False

    def attempt_set_variant_info_from_pre_existing_imported_allele_info(self) -> ImportedAlleleInfo:
        """
        Link to ImportedAlleleInfo (if haven't already), then update classification with any existing ImportedAlleleInfo
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
            :param user: The user creating this Classification
            :param lab: The lab the recod will be created under
            :param lab_record_id: The lab specific id for this record, will be auto-generated if not provided
            :param data: Initial data to populate with
            :param save: Should be True unless we're in test mode
            :param source: see patch_value
            :param make_fields_immutable: see patch_value
            :param populate_with_defaults: if True with populate data with defaults as defined by the EvidenceKeys overrides
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

        record.patch_value(data,
                           user=user,
                           clear_all_fields=True,
                           source=source,
                           save=save,
                           make_patch_fields_immutable=make_fields_immutable,
                           initial_data=True)
        if save:
            try:
                if making_new_id:
                    record.lab_record_id = 'vc' + str(record.id)
                    record.save()
            except:
                pass
                # in case of collision, leave the uuid

        return record

    def save(self, **kwargs):
        """ The post_save event is fired after save(), for other logic that depends on Classifications
            but doesn't work with the Classification object itself.
        """
        fix_permissions = kwargs.pop('fix_permissions', None) or self.id is None
        super().save(**kwargs)

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
        # Do a case insensitive check for each value against the key and any aliases
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
                if settings.VARIANT_CLASSIFICATION_REQUIRE_OVERWRITE_NOTE and \
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
        return EvidenceKeyMap.instance(lab=self.lab)

    def process_entry(self, cell: VCDataCell, source: str):
        """
        Normalises the value and adds any validation
        """
        # parse all notes looking for references
        value = cell.value
        e_key = cell.e_key
        note = cell.note

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
            # only provide mandatory validation for internal labs
            # don't want to dirty up imported read only classifications with errors
            cell.add_validation(code=ValidationCode.MANDATORY, severity='error', message='Missing mandatory value')

        # if we got an array of values
        if isinstance(value, list):
            if not value:
                # if the array is empty, treat that as a null value
                value = None
            elif e_key.value_type != EvidenceKeyValueType.MULTISELECT and len(value) == 1:
                # if we're not a multi-select key and we have an array with 1 value in it
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

            # handle selects with multiple values (comma sep string, array etc)
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
                    if not value.startswith("http"):  # this makes beautiful soap angry thinking we're asking it to go to the URL
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
                    message = f"Invalid date (expect yyyy-mm-dd): {ve}"
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
            if settings.VARIANT_CLASSIFICATION_AUTOFUZZ_AGE:
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
            # in case we have multiple last published records, this will  ensure we only hav 1
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
                    ignore_if_only_patching: Optional[Set[str]] = None) -> ClassificationPatchResponse:
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
            :param remove_api_immutable: If True, immutability level (under variantgrid) is removed from all fields. Requires source: SubissionSource.VariantGrid
            :param initial_data: if True, divides c.hgvs to
            :param revalidate_all: if True, runs validation over all fields we have, otherwise only the values being patched
            :param ignore_if_only_patching: if provided, if only these fields are different in the patch, don't both to patch anything
            :returns: A dict with "messages" (validation errors, warnings etc) and "modified" (fields that actually changed value)
        """
        source = source or SubmissionSource.API

        patch_response = ClassificationPatchResponse()

        is_admin_patch = source.can_edit(SubmissionSource.VARIANT_GRID)
        key_dict: EvidenceKeyMap = self.evidence_keys

        use_evidence = VCDataDict(copy.deepcopy(self.evidence),
                                  evidence_keys=self.evidence_keys)  # make a deep copy so we don't accidentally mutate the data
        patch = VCDataDict(data=EvidenceMixin.to_patch(patch),
                           evidence_keys=self.evidence_keys)  # the patch we're going to apply ontop of the evidence

        # make sure gene symbol is uppercase
        # need to do it here because it might get used in c.hgvs
        gene_symbol_cell = patch[SpecialEKeys.GENE_SYMBOL]
        if gene_symbol_cell.provided:
            if gene_symbol := gene_symbol_cell.value:
                if isinstance(gene_symbol, str) and gene_symbol != gene_symbol.upper():
                    gene_symbol_cell.value = gene_symbol.upper()

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

        # if submitting via API treat null as {value:None, explain:None, notes:None} for known keys
        # so we clear out any previous values but still retain immutability
        # for non existent keys just leave as null (if it started as null)
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
                    patch_response.append_warning(key=key, code="immutable", message=f"Cannot change immutable value for {key} from {existing.value} to blank")
                    patched.wipe(WipeMode.POP)  # reject entire change if attempting to change immutable value
                else:
                    patched.wipe(WipeMode.SET_NONE)
            else:
                patched.raw = cell.diff(dest=existing, ignore_if_omitted={'immutable'})
                if ('value' in patched or 'explain' in patched) and not source.can_edit(existing_immutability):
                    patch_response.append_warning(key=key, code="immutable", message=f"Cannot change immutable value or explain for {key} from {existing.value} to {patched.value}")
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
            clinical_significance_choice = self.calc_clinical_significance_choice()
            self.clinical_significance = clinical_significance_choice
            pending_modification.clinical_significance = clinical_significance_choice

            patch_response.append_warning(code="patched", message='Patched changed values for ' + ', '.join(sorted(diffs_only_patch.keys())))

            if save:
                self.save()
                patch_response.saved = True
                # have to save the modification after saving the record in case this
                # is the first record
                if pending_modification:
                    pending_modification.classification = self
                    pending_modification.save()
                    assign_perm(ClassificationModification.get_read_perm(), self.lab.group, pending_modification)

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

    def publish_latest(self, user: User, share_level=None, debug_timer: DebugTimer = DebugTimer.NullTimer) -> bool:
        if not share_level:
            share_level = self.share_level_enum
        latest_edited = self.last_edited_version
        if not latest_edited:
            raise ValueError(f'VC {self.id} does not have a last edited version')
        return latest_edited.publish(share_level=share_level, user=user, vc=self, debug_timer=debug_timer)

    @property
    def last_edited_version(self) -> 'ClassificationModification':
        return ClassificationModification.objects.filter(classification=self, is_last_edited=True).first()

    def latest_modification_for_user(self, user: User, exclude_withdrawn: bool = True) -> 'ClassificationModification':
        return ClassificationModification.latest_for_user(user, classification=self,
                                                          exclude_withdrawn=exclude_withdrawn).first()

    @property
    def last_published_version(self) -> 'ClassificationModification':
        return ClassificationModification.objects.filter(classification=self, is_last_published=True)\
            .select_related('classification', 'classification__lab', 'classification__lab__organization')\
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

    def as_json(self, params: ClassificationJsonParams) -> dict:
        current_user = params.current_user
        include_data = params.include_data
        version = params.version
        flatten = params.flatten
        include_lab_config = params.include_lab_config
        include_messages = params.include_messages
        strip_complicated = params.strip_complicated
        api_version = params.api_version

        if version and not isinstance(version, ClassificationModification):
            version = self.modification_at_timestamp(version)

        include_data_bool = include_data or flatten

        title = self.lab.group_name + '/' + self.lab_record_id
        # have version and last_edited as separate as we might be looking at an older version
        # but still want to know last edited
        latest_modification = None
        last_published_version = None

        if version:
            if version.is_last_published:
                last_published_version = version
            if version.is_last_edited:
                latest_modification = version

        if not latest_modification:
            latest_modification = self.last_edited_version
        if not last_published_version:
            last_published_version = self.last_published_version

        can_write = self.can_write(user=current_user)

        if self.share_level == ShareLevel.CURRENT_USER.key:
            last_published_version = None

        last_published_timestamp = None
        lowest_share_level = self.lowest_share_level(current_user)
        if last_published_version:
            last_published_timestamp = last_published_version.created.timestamp()

        clinical_context = None
        if self.clinical_context:
            clinical_context = self.clinical_context.name

        content = {
            'id': self.id,
            'lab_record_id': self.lab_record_id,
            'institution_name': self.lab.organization.name,
            'lab_id': self.lab.group_name,
            'org_name': self.lab.organization.shortest_name,
            'lab_name': self.lab.name,
            'title': title,
            'publish_level': self.share_level,
            'version_publish_level': version.share_level if version else self.share_level,
            'version_is_published': version.published if version else None,
            'is_last_published': version.is_last_published if version else (True if bool(self.id) else None),
            'published_version': last_published_timestamp,
            'can_write': can_write,
            'can_write_latest': can_write,
            'clinical_context': clinical_context,
            'data': self.get_visible_evidence(self.evidence, lowest_share_level),
            'resolved_condition': self.condition_resolution_dict,
            'withdrawn': self.withdrawn
        }
        if latest_modification:
            content['version'] = latest_modification.created.timestamp()

        if self.id:
            # only get flag collection if we've saved a record
            content['flag_collection'] = self.flag_collection_safe.pk

        if last_published_version and latest_modification:
            content['has_changes'] = last_published_version.id != latest_modification.id

        if include_messages:
            content['messages'] = Classification.validate_evidence(self.evidence)

        if can_write and latest_modification:
            content['last_edited'] = latest_modification.created.timestamp()

        if version is None or version.is_last_edited:
            # show the editable version
            pass
        else:
            version_timestamp = version.created.timestamp()
            versioned_data = version.evidence
            content['title'] = content['title'] + '.' + str(version_timestamp)
            content['version'] = version_timestamp
            content['data'] = self.get_visible_evidence(versioned_data, lowest_share_level)
            if include_messages:
                content['messages'] = Classification.validate_evidence(versioned_data)
            content['can_write'] = False

        if include_lab_config:
            content['config'] = EvidenceKey.merge_config(self.lab.classification_config,
                                                         self.lab.organization.classification_config)

        if include_messages:
            content['messages'] = (content['messages'] or []) + self.current_state_validation(
                user=current_user).flat_messages

        # Summary
        data = content['data']

        content["allele"] = self.get_allele_info_dict()

        if self.sample:
            content["sample_id"] = self.sample.pk
            content["sample_name"] = self.sample.name

        if include_data is not None and not isinstance(include_data, bool):
            if isinstance(include_data, str):
                include_data = [include_data]
            data_copy = {}
            for key in include_data:
                if key in data:
                    data_copy[key] = data[key]
                else:
                    data_copy[key] = None
            data = data_copy
        else:
            # do a copy just in case, want to make sure
            # that there are no errant saves going on
            data = content['data'].copy()

        if data and params.fix_data_types:
            # currently just fixes data types to multiselect array
            # should fix more data types and move the logic out of as_json so it can be done for other datatypes
            e_keys = EvidenceKeyMap.instance()
            for key, blob in data.items():
                if e_keys.get(key).value_type == EvidenceKeyValueType.MULTISELECT:
                    value = blob.get('value')
                    if isinstance(value, str):
                        blob['value'] = [s.strip() for s in value.split(',')]

        for key, blob in data.items():
            if strip_complicated:
                remove_keys = set(blob.keys()) - {'value', 'note', 'explain', 'db_refs'}
                for remove_key in remove_keys:
                    blob.pop(remove_key)
            elif not include_messages:
                blob.pop('validation', None)

        if flatten:
            data = Classification.flatten(data, ignore_none_values=True)

        if include_data_bool:
            content['data'] = data
        else:
            content.pop('data', None)

        if api_version >= 2:
            META_KEYS = [
                'id', 'lab_record_id', 'institution_name', 'lab_id', 'lab_name', 'title',
                'user', 'published_version', 'can_write', 'can_write_latest',
                'clinical_context', 'flag_collection', 'has_changes', 'version',
                'last_edited'
            ]
            meta = {}
            for meta_key in META_KEYS:
                if meta_key in content:
                    meta[meta_key] = content.pop(meta_key)
            if 'publish_level' in content:
                content['publish'] = content.pop('publish_level')

            content['meta'] = meta
            content['id'] = self.lab.group_name + '/' + self.lab_record_id

        if params.hardcode_extra_data:
            for key, value in params.hardcode_extra_data.items():
                content[key] = value

        return content

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
                value = v
            else:
                value = {'value': "(hidden)", 'hidden': True}
            visible_evidence[k] = value
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
                    (preferred_build := allele_info[self.get_genome_build()]) and \
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
    def get_classifications_qs(user: User, clinical_significance_list: Iterable[str] = None) -> QuerySet:
        vcm_qs = ClassificationModification.latest_for_user(user, published=True)
        if clinical_significance_list:
            vcm_qs = vcm_qs.filter(clinical_significance__in=clinical_significance_list)
        return Classification.objects.filter(pk__in=vcm_qs.values('classification'))

    @staticmethod
    def get_variant_q(user: User, genome_build: GenomeBuild, clinical_significance_list: Iterable[str] = None) -> Q:
        """ returns a Q object filtering variants to those with a PUBLISHED classification
            (optionally classification clinical significance in clinical_significance_list """
        vc_qs = Classification.get_classifications_qs(user, clinical_significance_list)
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
            variant_qs = Variant.objects.filter(variantallele__allele=OuterRef("classification__variant__variantallele__allele"),
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
            va_set = self.variant.allele.variantallele_set
            variant_allele: VariantAllele = va_set.filter(genome_build=genome_build).first()
            if variant_allele:
                return variant_allele.variant
        return None

    def get_variant_annotation(self, variant_annotation_version: VariantAnnotationVersion) -> Optional[Variant]:
        variant_annotation = None
        variant = self.get_variant_for_build(variant_annotation_version.genome_build)
        if variant:
            variant_annotation = variant.variantannotation_set.filter(version=variant_annotation_version).first()
        return variant_annotation

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

    def _generate_c_hgvs_extra(self, genome_build: GenomeBuild) -> HGVSNameExtra:
        variant = self.get_variant_for_build(genome_build)
        hgvs_matcher = HGVSMatcher(genome_build=genome_build)

        c_hgvs: HGVSNameExtra = HGVSNameExtra()
        if variant:
            transcript_id = None
            try:
                transcript_id = self.transcript
                c_hgvs = hgvs_matcher.variant_to_c_hgvs_extra(variant, transcript_id)
            except Exception:
                # can't map between builds
                report_exc_info(extra_data={
                    "genome_build": genome_build.name,
                    "variant": str(variant),
                    "transcript_id": transcript_id
                })
        return c_hgvs

    def get_c_hgvs(self, genome_build: GenomeBuild, use_full: bool = False) -> Optional[str]:
        if genome_build == genome_build.grch37():
            return self.chgvs_grch37 if not use_full else self.chgvs_grch37_full
        if genome_build == genome_build.grch38():
            return self.chgvs_grch38 if not use_full else self.chgvs_grch38_full
        return self._generate_c_hgvs_extra(genome_build).format()

    def __str__(self) -> str:
        genome_build = GenomeBuildManager.get_current_genome_build()
        cached_c_hgvs = self.get_c_hgvs(genome_build=genome_build)
        if not cached_c_hgvs:
            cached_c_hgvs = self.get(SpecialEKeys.C_HGVS)

        clinical_significance = self.get_clinical_significance_display() or "Unclassified"
        return f"({str(self.id)}) {cached_c_hgvs} {clinical_significance}"

    @staticmethod
    def check_can_create_no_classification_via_web_form(_user: User):
        if not settings.VARIANT_CLASSIFICATION_WEB_FORM_CREATE_ALLOW_NO_VARIANT:
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
    # it's useful if we have to rework permissions though
    share_level = models.CharField(max_length=16, choices=ShareLevel.choices(), null=True, blank=True)

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

    @staticmethod
    def column_name_for_build(genome_build: GenomeBuild, suffix: str = 'c_hgvs'):

        build_str: str
        if genome_build.is_equivalent(GenomeBuild.grch37()):
            build_str = 'grch37'

        elif genome_build.is_equivalent(GenomeBuild.grch38()):
            build_str = 'grch38'
        else:
            raise ValueError(f'No cached column for genome build {genome_build.pk}')
        return f'classification__allele_info__{build_str}__{suffix}'

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
        :param kwargs: Any other parameters will be used for filtering
        :return:
        """
        qs = get_objects_for_user(user, cls.get_read_perm(), klass=cls, accept_global_perms=True)
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
    def publish(self, share_level: Union[ShareLevel, str, int], user: User, vc: 'Classification', debug_timer: DebugTimer = DebugTimer.NullTimer) -> bool:
        """
        :param share_level: The share level we want to publish as
        :param user: The user who initiated the publishing
        :param vc: The variant classification we're publishing - even though it's redundant but we get the same object
        instance and don't reload it from the database
        :param debug_timer: for timing the event
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

        vc.save()

        debug_timer.tick("Published Modification")

        classification_post_publish_signal.send(
            sender=Classification,
            classification=vc,
            previously_published=previously_published,
            previous_share_level=old_share_level,
            newly_published=self,
            user=user,
            debug_timer=debug_timer)

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
        Determines as far as a user is concerned if two versions of a classification are equal based on evidence key data and share level.
        Would typically be used to determine if should show the column in history.
        Raises an error if run on different classification ids.
        :param other: Another modification
        :param care_about_explains: True if a change in only the explain should cause records not to be considered the same
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


class ClassificationConsensus:

    def __init__(self, variant: Variant, user: User, copy_from: Optional[ClassificationModification] = None):
        variants = [variant]
        allele = variant.allele
        self.user = user
        if copy_from:
            self.vcm = copy_from
            copy_from.check_can_view(user)
        else:
            self.vcm: Optional[ClassificationModification] = None
            if allele:
                variants = allele.variants
            vcms = list(ClassificationModification.latest_for_user(user=user, variant=variants, published=True).filter(classification__lab__external=False).all())
            if vcms:
                vcms.sort(key=lambda vcm: vcm.curated_date_check)
                self.vcm = vcms[-1]

    @property
    def is_valid(self) -> bool:
        return self.vcm is not None

    @cached_property
    def consensus_patch(self) -> VCPatch:
        keys = EvidenceKeyMap.instance()
        if not self.vcm:
            return {}

        evidence = self.vcm.published_evidence
        consensus: Dict[str, Any] = {}
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


@dataclass
class ClassificationDateType:
    date: datetime
    name: Optional[str] = None

    @property
    def is_fallback(self) -> bool:
        return self.name == "Created"


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

    def convert_date(self, evidence_key):
        if date_str := self._modification.get(evidence_key):
            try:
                return Classification.to_date(date_str).astimezone(self.timezone)
            except:
                pass
        return None

    @cached_property
    def curated_date(self) -> Optional[datetime]:
        return self.convert_date(SpecialEKeys.CURATION_DATE)

    @cached_property
    def sample_date(self) -> Optional[datetime]:
        return self.convert_date(SpecialEKeys.SAMPLE_DATE)

    @cached_property
    def curated_verified_date(self) -> Optional[datetime]:
        return self.convert_date(SpecialEKeys.CURATION_VERIFIED_DATE)

    @cached_property
    def created_date(self) -> datetime:
        return self._modification.classification.created.astimezone(self.timezone)

    @cached_property
    def epoch(self):
        return datetime.utcfromtimestamp(0).astimezone(self.timezone)

    @cached_property
    def relevant_date(self) -> ClassificationDateType:
        """
        Return the most relevant date (for reporting) along with a name for the date if it isn't the Sample Date
        """
        if date := self.curated_date:
            return ClassificationDateType(date)  # don't provide label for curated as it's the default
        elif date := self.curated_verified_date:
            return ClassificationDateType(date, "Verified")
        elif date := self.sample_date:
            return ClassificationDateType(date, "Sample Date")
        return ClassificationDateType(self.created_date, "Created")

    @property
    def name(self) -> str:
        return self.relevant_date.name

    @property
    def date(self) -> datetime:
        return self.relevant_date.date

    def __lt__(self, other: 'CuratedDate') -> bool:
        epoch = self.epoch

        my_relevant_date = self.relevant_date
        other_relevant_date = other.relevant_date

        # fallback dates (e.g. Created) are seen as before all other dates
        # can almost think of it as null... but not 100% null
        if my_relevant_date.is_fallback != other_relevant_date.is_fallback:
            return my_relevant_date.is_fallback

        # both dates aren't fallback (or both are) if they're diff return the earliest
        if my_relevant_date.date != other_relevant_date.date:
            return my_relevant_date.date < other_relevant_date.date

        # Most relevant dates are the same, now go down the chain until we find a difference

        def direction(date_1, date_2):
            nonlocal epoch
            date_1 = date_1 or epoch
            date_2 = date_2 or epoch

            if date_1 == date_2:
                return 0
            elif date_1 < date_2:
                return -1
            else:
                return 1

        if diff := direction(self.curated_date, other.curated_date):
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
