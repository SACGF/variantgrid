from dataclasses import dataclass
from functools import cached_property
from typing import Optional

from auditlog.context import set_extra_data
from django.contrib.auth.models import User
from django.dispatch import receiver
from django.urls import reverse

from classification.enums import ClassificationResultValue, ShareLevel, SpecialEKeys, TriageDecision
from classification.enums.overlaps_enums import TriageStatus
from classification.models import (
    Classification,
    ClassificationFlagTypes,
    ClassificationGrouping,
    ClassificationModification,
    ClassificationSummaryCacheObj,
    EvidenceKey,
    EvidenceKeyMap,
    OverlapContribution,
    classification_flag_types,
    classification_grouping_summary_signal,
    classification_post_publish_signal,
)
from classification.services.overlaps_services import OverlapServices
from eventlog.notifications import NotificationBuilder
from library.django_utils import get_url_from_view_path
from library.guardian_utils import admin_bot
from library.utils import get_timer
from snpdb.enums import ClinicalSignificance


@receiver(classification_post_publish_signal, sender=Classification)
def turn_off_unsubmitted_edits(sender, classification, previously_published, newly_published, user, **kwargs):  # pylint: disable=unused-argument
    """
    Removes outstanding edit flags if there aren't any left.
    This is a bit messy as regular code inserts it, but a signal closes it.
    (somewhat redundant to Classification.insert() but if publish is called without an insert, the outstanding edits flag wont be closed)
    """
    classification.flag_collection_safe.ensure_resolution(
        flag_type=classification_flag_types.classification_outstanding_edits,
        resolution='closed'
    )


@receiver(classification_grouping_summary_signal, sender=ClassificationGrouping)
def summary_changed(sender, instance: ClassificationGrouping, old_summary: ClassificationSummaryCacheObj, new_summary: ClassificationSummaryCacheObj, **kwargs):

    for value_type in (ClassificationResultValue.ONC_PATH, ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE):
        if old_summary.value_for_value_type(value_type) != new_summary.value_for_value_type(value_type):
            if oc := OverlapContribution.objects.filter(classification_grouping=instance, value_type=value_type).first():
                if oc.triage_state_obj.status == TriageStatus.REVIEWED_WILL_FIX:
                    # Note any value change will cause AMENDED, even if it wasn't on the agreed upon value
                    with set_extra_data({
                        "actor": admin_bot(),
                        "remote_addr": None,
                        "remote_port": None,
                    }):
                        oc.triage_state_obj = TriageDecision(status=TriageStatus.AMENDED)
                        oc.save()

                    for overlap in oc.overlaps:
                        OverlapServices.recalc_overlap(overlap)


# THE BELOW ARE FOR SINGLE CLASSIFICATION FLAGS
# these flags will be removed in future, either replaced with a notification on the grouping that the value has changed
# or with nothing at all

@dataclass
class _SignificanceChange:
    previously_published: ClassificationModification
    newly_published: ClassificationModification
    evidence_key: EvidenceKey

    @cached_property
    def previous_value(self) -> Optional[str]:
        return self.previously_published.get(self.evidence_key.key)

    @cached_property
    def new_value(self) -> Optional[str]:
        return self.newly_published.get(self.evidence_key.key)

    def __bool__(self):
        return self.is_significant

    @cached_property
    def is_significant(self) -> bool:
        return ClinicalSignificance.is_significant_change(self.previous_value, self.new_value)

    @cached_property
    def previous_label(self) -> str:
        return self.evidence_key.pretty_value(self.previous_value) or "No Data"

    @cached_property
    def new_label(self) -> str:
        return self.evidence_key.pretty_value(self.new_value) or "No Data"

    @property
    def comment(self):
        return f'{self.evidence_key.pretty_label} changed from {self.previous_label} to {self.new_label}'

    @property
    def pending_close_comment(self):
        return f'{self.evidence_key.pretty_label} changed to {self.new_label}'

    def raise_flag(self, classification: Classification):
        classification.flag_collection_safe.add_flag(
            flag_type=classification_flag_types.significance_change,
            comment=self.comment,
            data={
                "e_key": self.evidence_key.key,
                "old_value": self.previous_value,
                "new_value": self.new_value
            }
        )


@receiver(classification_post_publish_signal, sender=Classification)
def clinical_significance_change_check(
        sender,
        classification:Classification,
        previously_published: ClassificationModification,
        newly_published: ClassificationModification,
        user: User,
        **kwargs):  # pylint: disable=unused-argument
    """
    Closes unshared flag
    Used to handle pending changes, but that is now done with triage states isntead of flags on individual classifications
    """

    # TODO move this into classification_hooks_share_flags
    if classification.share_level_enum.index > ShareLevel.ORGANISATION.index:
        classification.flag_collection_safe.close_open_flags_of_type(
            flag_type=classification_flag_types.unshared_flag
        )