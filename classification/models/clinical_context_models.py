from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Iterable, Set

import django.dispatch
from django.conf import settings
from django.contrib.auth.models import User
from django.db import models, transaction
from django.db.models.deletion import CASCADE
from django.db.models.query import QuerySet
from django.dispatch.dispatcher import receiver
from django.utils.timezone import now
from django_extensions.db.models import TimeStampedModel
from lazy import lazy

from classification.enums import ShareLevel, SpecialEKeys
from classification.enums.clinical_context_enums import ClinicalContextStatus
from classification.models.classification import Classification, \
    ClassificationModification
from classification.models.classification_import_run import ClassificationImportRun, \
    classification_imports_complete_signal
from classification.models.flag_types import classification_flag_types
from flags.models.models import FlagsMixin, FlagCollection, FlagTypeContext, \
    flag_collection_extra_info_signal, FlagInfos
from library.django_utils import get_url_from_view_path
from library.log_utils import report_message, NotificationBuilder
from snpdb.models import Variant, Lab
from snpdb.models.models_variant import Allele

clinical_context_signal = django.dispatch.Signal(providing_args=["clinical_context", "status", "is_significance_change", "cause"])

CS_TO_NUMBER = {
    'B': 1,
    'LB': 1,
    'VUS': 2,
    'VUS_A': 2,
    'VUS_B': 2,
    'VUS_C': 2,
    'LP': 3,
    'P': 3,
    'R': 4
}
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
        if self in (DiscordanceLevel.CONCORDANT_AGREEMENT, DiscordanceLevel.CONCORDANT_DIFF_VUS):
            return "overlap-agreement"
        if self == DiscordanceLevel.CONCORDANT_CONFIDENCE:
            return "overlap-confidence"
        if self == DiscordanceLevel.DISCORDANT:
            return "overlap-discordant"
        return "overlap-single"


@dataclass
class DiscordanceStatus:
    level: DiscordanceLevel
    lab_count: int
    lab_count_all: int
    counted_classifications: int

    @staticmethod
    def calculate(modifications: Iterable[ClassificationModification]) -> 'DiscordanceStatus':
        cs_scores: Set[int] = set()  # clin sig to score, all VUSs are 3
        cs_vuses: Set[int] = set()  # all the different VUS A,B,C values
        cs_values: Set[str] = set()   # all the different clinical sig values
        shared_labs: Set[Lab] = set()   # all the sharing labs
        all_labs: Set[Lab] = set()   # all labs involved (sharing and not sharing, only sharing records determine the status)

        level: Optional[DiscordanceLevel]
        counted_classifications = 0
        for vcm in modifications:
            all_labs.add(vcm.classification.lab)
            clin_sig = vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            if vcm.share_level_enum.is_discordant_level and vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE) is not None:
                if strength := CS_TO_NUMBER.get(clin_sig):
                    counted_classifications += 1
                    shared_labs.add(vcm.classification.lab)
                    cs_scores.add(strength)
                    cs_values.add(clin_sig)
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

        return DiscordanceStatus(level=level, lab_count=len(shared_labs), lab_count_all=len(all_labs), counted_classifications=counted_classifications)


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

    def calculate_status(self) -> str:
        """
        deprecated, try to use DiscordanceLevels instead
        """
        status = self.calculate_discordance_status()
        level = status.level
        if level == DiscordanceLevel.DISCORDANT:
            return ClinicalContextStatus.DISCORDANT
        else:
            return ClinicalContextStatus.CONCORDANT

    def calculate_discordance_status(self) -> DiscordanceStatus:
        return DiscordanceStatus.calculate(self.classification_modifications)

    @lazy
    def relevant_classification_count(self) -> int:
        return len([vcm for vcm in self.classification_modifications if CS_TO_NUMBER.get(vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))])

    def is_discordant(self):
        return self.status == ClinicalContextStatus.DISCORDANT

    @lazy
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
        :param cause_code: A finite list of codes as to what triggered the discordance (not currently used)
        """
        old_status = self.status
        new_status = self.calculate_status()
        is_significance_change = old_status != new_status

        ongoing_import = ClassificationImportRun.ongoing_import()

        self.last_evaluation = {
            "date": now().timestamp(),
            "trigger": cause,
            "considers_vcms": [vcm.id_str for vcm in self.classification_modifications],
            "old_status": old_status,
            "new_status": new_status,
            "delayed": ongoing_import
        }

        cause = self.pending_cause or cause

        if ongoing_import and self.status != new_status:
            if new_status != self.pending_status:
                # only update the cause if we're going to end up with a new status
                self.pending_cause = cause
                self.pending_status = new_status

                allele_url = get_url_from_view_path(self.allele.get_absolute_url())
                nb = NotificationBuilder("ClinicalContext changed (delayed due to ongoing import)")
                nb.add_markdown(
                    f"ClinicalGrouping for allele <{allele_url}|{allele_url}> would change from {old_status} -> {new_status} but marked as pending due to {ongoing_import}")
                nb.send()

        else:
            # if doing the recalc live OR if delayed but the status has reverted to what it used to be
            # wipe out the old values
            self.status = new_status
            if self.pending_status:
                allele_url = get_url_from_view_path(self.allele.get_absolute_url())
                nb = NotificationBuilder("ClinicalContext changed-back (delayed due to ongoing import)")
                nb.add_markdown(
                    f"ClinicalGrouping for allele <{allele_url}|{allele_url}> changed back from {self.pending_status} -> {new_status} within {ongoing_import}, no notifications sent")
                nb.send()

            self.pending_cause = None
            self.pending_status = None

        self.save()

        if ongoing_import:
            # don't send out any signals if calculation is delayed, we don't want to trigger anything until we complete
            return

        clinical_context_signal.send(sender=ClinicalContext, clinical_context=self, status=new_status, is_significance_change=is_significance_change, cause=cause)
        # might be worth checking to see if the clinical context has changed against post signal

        if settings.DISCORDANCE_ENABLED:
            if new_status == ClinicalContextStatus.DISCORDANT:
                self.flag_collection_safe.get_or_create_open_flag_of_type(
                    flag_type=classification_flag_types.clinical_context_discordance,
                    reopen=True
                )

            else:
                self.flag_collection_safe.close_open_flags_of_type(
                    flag_type=classification_flag_types.clinical_context_discordance
                )

            if is_significance_change and (old_status or new_status == ClinicalContextStatus.DISCORDANT):
                report_message('Allele ID clinical grouping change', extra_data={
                    'target': f'Allele ID {self.allele_id} {old_status} -> {new_status}'
                })

        else:
            self.flag_collection_safe.close_open_flags_of_type(
                flag_type=classification_flag_types.clinical_context_discordance,
                comment='Discordance functionality has been disabled'
            )

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
    def classifications_qs(self) -> QuerySet:
        return Classification.objects.filter(clinical_context=self, share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS, withdrawn=False)

    @property
    def classifications_associated_qs(self) -> QuerySet:
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
