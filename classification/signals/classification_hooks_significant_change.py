from dataclasses import dataclass
from functools import cached_property
from typing import Optional

from django.dispatch import receiver
from django.urls import reverse
from django_messages.admin import User

from classification.enums import ShareLevel, SpecialEKeys, ClinicalSignificance
from classification.models import classification_post_publish_signal, Classification, classification_flag_types, \
    EvidenceKey, ClassificationFlagTypes, ClassificationModification, EvidenceKeyMap
from library.django_utils import get_url_from_view_path
from library.log_utils import NotificationBuilder
from library.utils import get_timer


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
    Raises a clinical significance change flag if the clinical significance has changed after publishing
    """

    # TODO move this into classification_hooks_share_flags
    if classification.share_level_enum.index > ShareLevel.INSTITUTION.index:
        classification.flag_collection_safe.close_open_flags_of_type(
            flag_type=classification_flag_types.unshared_flag
        )

    if previously_published and previously_published.share_level in ShareLevel.DISCORDANT_LEVEL_KEYS:

        if somatic_clin_sig_change := _SignificanceChange(
            previously_published=previously_published,
            newly_published=newly_published,
            evidence_key=EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)
        ):
            somatic_clin_sig_change.raise_flag(classification)

        if classification_change := _SignificanceChange(
            previously_published=previously_published,
            newly_published=newly_published,
            evidence_key=EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        ):
            classification_change.raise_flag(classification)
            close_message = f"Classification was changed to {classification_change.new_label}"

            if pending_change_flag := classification.flag_collection_safe.get_flag_of_type(classification_flag_types.classification_pending_changes):
                if pending_change_value := pending_change_flag.data.get(ClassificationFlagTypes.CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY):

                    class_url = get_url_from_view_path(
                        reverse('view_classification', kwargs={'classification_id': classification.id}),
                    )

                    pending_change_label = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(pending_change_value)
                    if classification_change.new_label != pending_change_value:
                        close_message += f", expected {pending_change_label}"

                    # is_agreed_change = newly_published != pending_change_value
                    nb = NotificationBuilder("Pending Change")
                    nb.add_markdown(
                        f"Classification <{class_url}|{classification.friendly_label}> had a pending change and has now been updated:\n"
                        f"*{classification_change.previous_label}* - old clinical significance\n"
                        f"*{pending_change_label}* - pending change\n"
                        f"*{classification_change.new_label}* - new clinical significance")
                    nb.send()

            # note we don't care if the clin sig was changed to agreed upon value per the flag
            # just that it was changed
            classification.flag_collection_safe.close_open_flags_of_type(
                flag_type=classification_flag_types.classification_pending_changes,
                comment=close_message
            )

    get_timer().tick("Update share flags")
