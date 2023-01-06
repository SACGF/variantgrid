from dataclasses import dataclass
from enum import Enum
from functools import cached_property
from typing import List, Optional, Iterable, Set, Dict

import django.dispatch
from django.conf import settings
from django.contrib.auth.models import User
from django.db import models, transaction
from django.db.models.deletion import CASCADE
from django.db.models.query import QuerySet
from django.dispatch.dispatcher import receiver
from django.utils.timezone import now
from django_extensions.db.models import TimeStampedModel

from classification.enums import ShareLevel, SpecialEKeys
from classification.enums.clinical_context_enums import ClinicalContextStatus
from classification.models.classification import Classification, \
    ClassificationModification
from classification.models.classification_import_run import ClassificationImportRun, \
    classification_imports_complete_signal
from flags.models import Flag, FlagStatus
from flags.models.models import FlagsMixin, FlagCollection, FlagTypeContext, \
    flag_collection_extra_info_signal, FlagInfos
from library.django_utils import get_url_from_view_path
from library.log_utils import NotificationBuilder
from library.utils import invalidate_cached_property
from snpdb.models import Variant, Lab
from snpdb.models.models_variant import Allele

clinical_context_signal = django.dispatch.Signal()  # args: "clinical_context", "status", "is_significance_change", "cause"

# TODO, consider moving this into the clinical significance evidence key options rather than hardcoded
SPECIAL_VUS = {
    'VUS_A': 1,
    'VUS_B': 2,
    'VUS_C': 3
}


class DiscordanceLevel(str, Enum):
    """
    Values are assigned based on how big the differences are considered when it comes to discordance
    e.g. difference between Benign and Likely Benign isn't as important as the difference between Likely Benign and VUS etc
    """

    NO_ENTRIES = 'no_entries'
    SINGLE_SUBMISSION = 'single_submission'  # single shared submission
    CONCORDANT_AGREEMENT = 'concordant_agreement'  # complete agreement
    CONCORDANT_DIFF_VUS = 'concordant_agreement_diff_vus'  # VUS-A, VUS-B vs VUS-C
    CONCORDANT_CONFIDENCE = 'concordant_confidence'  # Benign vs Likely Benign
    DISCORDANT = 'discordant'

    @property
    def label(self) -> str:
        if self == DiscordanceLevel.CONCORDANT_AGREEMENT:
            return "Concordant (Agreement)"
        if self == DiscordanceLevel.CONCORDANT_DIFF_VUS:
            return "Concordant (Agreement Differing VUS)"
        if self == DiscordanceLevel.CONCORDANT_CONFIDENCE:
            return "Concordant (Confidence)"
        if self == DiscordanceLevel.NO_ENTRIES:
            return "No Shared Submissions"
        if self == DiscordanceLevel.SINGLE_SUBMISSION:
            return "Single Shared Submission"
        if self == DiscordanceLevel.DISCORDANT:
            return "Discordant"
        return "Unknown"

    @property
    def css_class(self) -> str:
        return f"overlap-{self.value}"


@dataclass(frozen=True)
class _DiscordanceCalculationRow:
    lab: Lab
    clinical_significance: str
    shared: bool


@dataclass
class DiscordanceStatus:
    """
    TODO rename DiscordanceStatus to OverlapSummary

    There's a bit of confusion due to evolution of discordance status
    DiscordanceStatus: gives absolute full context about a clinical grouping
    DiscordanceLevel: an enum that covers all states, e.g. no submissions, single submission, overlap confidence
    DiscordanceReportStatus: how the user left the DiscordanceReport
    ClinicalContextStatus: saved on the clinical grouping and saved to the database, maybe it should be repalced with DiscordanceLevel?
    """
    level: DiscordanceLevel
    lab_count: int
    lab_count_all: int
    counted_classifications: int
    pending_concordance: bool = False
    discordance_report: Optional['DiscordanceReport'] = None
    has_ignored_clin_sigs: bool = False

    def __str__(self):
        if self.pending_concordance:
            return "Pending Concordance"
        elif self.discordance_report and self.discordance_report.resolution == 'D':
            return "Continued Discordance"
        elif self.discordance_report and self.discordance_report.resolution is None:
            return "Active Discordance"
        else:
            return self.level.label

    @property
    def sort_order(self):
        if self.pending_concordance:
            return 5
        elif self.discordance_report and self.discordance_report.resolution == 'D':
            return 6
        elif self.discordance_report and self.discordance_report.resolution is None:
            return 7
        elif self.level == DiscordanceLevel.CONCORDANT_CONFIDENCE:
            return 4
        elif self.level == DiscordanceLevel.CONCORDANT_DIFF_VUS:
            return 3
        elif self.level == DiscordanceLevel.CONCORDANT_AGREEMENT:
            return 2
        else:
            # single submissions, no shred submissions, shouldn't appear in overlaps page
            # and level discordance has already been taken care of
            return 1

    def __lt__(self, other):
        return self.sort_order < other.sort_order

    @property
    def css_class(self):
        if self.pending_concordance:
            return 'overlap-pending_concordance'
        else:
            return self.level.css_class

    @property
    def is_discordant(self):
        return self.level == DiscordanceLevel.DISCORDANT

    @staticmethod
    def cs_buckets():
        from classification.models import EvidenceKeyMap
        return EvidenceKeyMap.clinical_significance_to_bucket()

    @staticmethod
    def calculate(modifications: Iterable[ClassificationModification], discordance_report: Optional['DiscordanceReport'] = None) -> 'DiscordanceStatus':
        base_status = DiscordanceStatus._calculate_rows([
            _DiscordanceCalculationRow(
                lab=vcm.classification.lab,
                clinical_significance=vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE),
                shared=vcm.share_level_enum.is_discordant_level
            ) for vcm in modifications
        ])
        if not base_status.is_discordant:
            # no need to check flags for pending changes if we're not discordant when looking at all entries
            # the common case by far
            return base_status

        # we're discordant, check to see if there are pending changes AND if those pending changes would bring us into discordance
        # if so, just mark "pending_concordance" as true
        from classification.models import ClassificationFlagTypes, classification_flag_types
        flag_collections = {mod.classification.flag_collection_id: mod for mod in modifications}
        if pending_flags := Flag.objects.filter(collection__in=flag_collections.keys(),
                                                flag_type=classification_flag_types.classification_pending_changes,
                                                resolution__status=FlagStatus.OPEN):

            clin_sig_overrides: Dict[int, str] = {}
            for flag in pending_flags:
                clin_sig_overrides[flag_collections[flag.collection_id].pk] = flag.data.get(
                    ClassificationFlagTypes.CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY) or None

            override_status = DiscordanceStatus._calculate_rows([
                _DiscordanceCalculationRow(
                    lab=vcm.classification.lab,
                    clinical_significance=clin_sig_overrides.get(vcm.pk) if vcm.pk in clin_sig_overrides else vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE),
                    shared=vcm.share_level_enum.is_discordant_level
                ) for vcm in modifications
            ])
            if not override_status.is_discordant:
                base_status.pending_concordance = True

        base_status.discordance_report = discordance_report
        return base_status

    @staticmethod
    def _calculate_rows(rows: Iterable[_DiscordanceCalculationRow]) -> 'DiscordanceStatus':
        cs_scores: Set[int] = set()  # clin sig to score, all VUSs are 3
        cs_vuses: Set[int] = set()  # all the different VUS A,B,C values
        cs_values: Set[str] = set()   # all the different clinical sig values
        shared_labs: Set[Lab] = set()   # all the sharing labs
        all_labs: Set[Lab] = set()   # all labs involved (sharing and not sharing, only sharing records determine the status)
        ignored_clin_sigs: Set[str] = set()  # all clin sigs that we came across but

        level: Optional[DiscordanceLevel]
        counted_classifications = 0
        for row in rows:
            all_labs.add(row.lab)
            clin_sig = row.clinical_significance
            if clin_sig and row.shared:
                # if vcm.share_level_enum.is_discordant_level and vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE) is not None:
                if strength := DiscordanceStatus.cs_buckets().get(clin_sig):
                    counted_classifications += 1
                    shared_labs.add(row.lab)
                    cs_scores.add(strength)
                    cs_values.add(clin_sig)
                else:
                    shared_labs.add(row.lab)
                    ignored_clin_sigs.add(clin_sig)

                if vus_special := SPECIAL_VUS.get(clin_sig):
                    cs_vuses.add(vus_special)

        if counted_classifications == 0:
            level = DiscordanceLevel.NO_ENTRIES

        elif counted_classifications == 1:
            level = DiscordanceLevel.SINGLE_SUBMISSION

        elif len(cs_scores) > 1:
            level = DiscordanceLevel.DISCORDANT

            # not discordant if we've reached here, see if we have multiple VUS kinds
        elif len(cs_vuses) > 1:
            # importantly you can have a VUS vs VUS_A and still be in agreement
            # it's only if you have more than one of VUS_A,B,C that cs_vuses will have multiple values
            level = DiscordanceLevel.CONCORDANT_DIFF_VUS

        elif len(cs_values) > 1 and not (len(cs_scores) == 1 and 2 in cs_scores and len(cs_vuses) == 1):
            # check to make sure we don't have VUS vs VUS_A, or VUS vs VUS_B which should be considered agreemetn

            # importantly you can have a VUS vs VUS_A and still be in agreement
            # it's only if you have more than one of VUS_A,B,C that cs_vuses will have multiple values
            level = DiscordanceLevel.CONCORDANT_CONFIDENCE

        else:
            level = DiscordanceLevel.CONCORDANT_AGREEMENT

        return DiscordanceStatus(
            level=level,
            lab_count=len(shared_labs),
            lab_count_all=len(all_labs),
            counted_classifications=counted_classifications,
            has_ignored_clin_sigs=bool(ignored_clin_sigs)
        )


class ClinicalContextRecalcTrigger(Enum):
    ADMIN = "admin"
    VARIANT_SET = "variant_set"
    CLINICAL_GROUPING_SET = "clinical_grouping"
    WITHDRAW_DELETE = "withdraw_delete"
    SUBMISSION = "submission"
    DELAYED = "delayed"
    OTHER = "other"


class ClinicalContext(FlagsMixin, TimeStampedModel):
    default_name = 'default'

    # variant = models.ForeignKey(Variant, on_delete=CASCADE)
    allele = models.ForeignKey(Allele, on_delete=CASCADE)
    name = models.TextField(null=True, blank=True)
    status = models.TextField(choices=ClinicalContextStatus.CHOICES, null=True, blank=True)
    last_evaluation = models.JSONField(null=True, blank=True)

    pending_cause = models.TextField(null=True, blank=True)
    pending_status = models.TextField(choices=ClinicalContextStatus.CHOICES, null=True, blank=True)

    class Meta:
        unique_together = ('allele', 'name')

    def flag_type_context(self) -> FlagTypeContext:
        return FlagTypeContext.objects.get(pk='clinical_context')

    def calculate_status(self) -> ClinicalContextStatus:
        """
        deprecated, try to use DiscordanceLevels instead
        (but be careful there is an 'invalidate' calculate_status)
        """
        level = self.discordance_status.level
        if level == DiscordanceLevel.DISCORDANT:
            return ClinicalContextStatus.DISCORDANT
        else:
            return ClinicalContextStatus.CONCORDANT

    @property
    def is_pending_concordance(self):
        return self.discordance_status.pending_concordance

    def calculate_discordance_status(self) -> DiscordanceStatus:
        return DiscordanceStatus.calculate(self.classification_modifications, self.latest_report)

    @cached_property
    def discordance_status(self) -> DiscordanceStatus:
        # TODO this is awkwardly named as it calculates (at least the first time)
        # Need to distinguish between this and database cached is_discordant
        return DiscordanceStatus.calculate(self.classification_modifications, self.latest_report)

    @cached_property
    def relevant_classification_count(self) -> int:
        return len([vcm for vcm in self.classification_modifications if DiscordanceStatus.cs_buckets().get(vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))])

    def is_discordant(self):
        # WARNING: THIS IS APPLIES TO CACHED STATUS ONLY
        return self.status == ClinicalContextStatus.DISCORDANT

    @cached_property
    def latest_report(self) -> Optional['DiscordanceReport']:
        from classification.models import DiscordanceReport
        return DiscordanceReport.latest_report(self)

    @transaction.atomic
    def recalc_and_save(self,
                        cause: str,
                        cause_code: ClinicalContextRecalcTrigger = ClinicalContextRecalcTrigger.OTHER):
        """
        Updates this ClinicalContext with the new status and applies flags where appropriate.
        :param cause: A human readable string to be showed to the users as to what started or stopped a discordance (if relevant)
        typically will be "Lab X submitted a new variant"
        :param cause_code: A finite list of codes as to what triggered the discordance
        """
        invalidate_cached_property(self, 'discordance_status')

        old_status: ClinicalContextStatus = self.status
        new_status: ClinicalContextStatus = self.calculate_status()
        is_significance_change = old_status != new_status
        is_simple_new = old_status is None and new_status == ClinicalContextStatus.CONCORDANT
        allele_url = get_url_from_view_path(self.allele.get_absolute_url())

        ongoing_import = ClassificationImportRun.ongoing_imports()

        self.last_evaluation = {
            "date": now().timestamp(),
            "trigger": cause,
            "considers_vcms": [vcm.id_str for vcm in self.classification_modifications],
            "old_status": old_status,
            "new_status": new_status,
            "delayed": ongoing_import
        }

        cause = self.pending_cause or cause or "Unknown cause"

        if ongoing_import:
            # always set pending cause, to ensure calculation is re-triggered at the end
            # because even if the status doesn't change, a new lab might be added to an existing discordance report
            # so we still need to recaluclate

            self.pending_cause = self.pending_cause or cause
            # Important - this pending status is not to do with "Pending Changes"
            was_already_pending = self.pending_status != self.status
            pending_or_old_status = self.pending_status or self.status

            if (not is_simple_new) and (new_status != pending_or_old_status):
                # only update the cause if we're going to end up with a new status
                self.pending_status = new_status
                if settings.DISCORDANCE_ENABLED:
                    if new_status != old_status:
                        nb = NotificationBuilder("PENDING: ClinicalContext changed)")
                        nb.add_markdown(
                            f":clock1: ClinicalGrouping for allele <{allele_url}|{allele_url}> would change from {old_status} -> {new_status} but marked as pending due to {ongoing_import}")
                        nb.send()
                    else:
                        nb = NotificationBuilder("PENDING: ClinicalContext changed-back")
                        nb.add_markdown(
                            f":clock430: ClinicalGrouping for allele <{allele_url}|{allele_url}> changed back from {self.pending_status} -> {new_status} within {ongoing_import}, no notifications sent")
                        nb.send()

            self.save()
        else:
            # if doing the recalc live OR if delayed but the status has reverted to what it used to be
            # wipe out the old values
            if settings.DISCORDANCE_ENABLED and is_significance_change and not is_simple_new:
                nb = NotificationBuilder("ClinicalContext changed")
                nb.add_markdown(
                    f":fire_engine: ClinicalGrouping for allele <{allele_url}|{allele_url}> changed from {old_status} -> {new_status} "
                    f"\nLab notifications should follow")
                nb.send()

            self.status = new_status
            self.pending_cause = None
            self.pending_status = None

            self.save()

            # clinical_context_signal is now in charge of applying all relevant flags to clinical context and classifications
            clinical_context_signal.send(sender=ClinicalContext, clinical_context=self, status=new_status, is_significance_change=is_significance_change, cause=cause)

    @property
    def is_default(self) -> bool:
        return self.name == ClinicalContext.default_name

    @staticmethod
    def default_group_for(variant: Variant) -> Optional['ClinicalContext']:
        if variant.allele:
            dg, _ = ClinicalContext.objects.get_or_create(allele=variant.allele, name=ClinicalContext.default_name)
            return dg
        print("Warning, classification variant doesn't have allele")
        return None

    @staticmethod
    def for_variant(variant: Variant) -> QuerySet:
        ClinicalContext.default_group_for(variant)  # ensure we always have the default group
        return ClinicalContext.objects.filter(allele=variant.allele).order_by('created')

    @property
    def classifications_qs(self) -> QuerySet[Classification]:
        return Classification.objects.filter(clinical_context=self, share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS, withdrawn=False)

    @property
    def classifications_associated_qs(self) -> QuerySet[Classification]:
        return Classification.objects.filter(clinical_context=self)

    @property
    def classification_modifications(self) -> List[ClassificationModification]:
        vcms = [vcm for vcm in (vc.last_published_version for vc in self.classifications_qs) if vcm is not None]
        return vcms

    def to_json(self) -> dict:
        return {
            'id': self.id,
            'name': self.name,
            'is_default': self.is_default
        }

    def __str__(self) -> str:
        if self.is_default:
            return 'Default Grouping for Allele'
        return self.name


@receiver(flag_collection_extra_info_signal, sender=FlagCollection)
def get_extra_info(flag_infos: FlagInfos, user: User, **kwargs):  # pylint: disable=unused-argument
    ccs = ClinicalContext.objects.filter(flag_collection__in=flag_infos.ids).select_related('allele')
    for cc in ccs:
        flag_infos.set_extra_info(cc.flag_collection_id, {
            'label': f'{cc.allele} - "{cc.name}"',
            'cc_id': cc.id
        }, source_object=cc)


@receiver(classification_imports_complete_signal, sender=ClassificationImportRun)
def import_complete(**kwargs):
    # this is called when there are no ongoing imports, find all the delayed clinical contexts, and calculate them
    for cc in ClinicalContext.objects.filter(pending_cause__isnull=False):
        # cause should automatically be loaded from pending cause anyway
        cc.recalc_and_save(cause=cc.pending_cause, cause_code=ClinicalContextRecalcTrigger.DELAYED)
