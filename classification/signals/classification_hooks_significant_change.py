from django.dispatch import receiver
from django.urls import reverse

from classification.enums import ShareLevel, SpecialEKeys, ClinicalSignificance
from classification.models import classification_post_publish_signal, Classification, classification_flag_types, \
    EvidenceKey, ClassificationFlagTypes
from library.django_utils import get_url_from_view_path
from library.log_utils import NotificationBuilder
from library.utils import DebugTimer


@receiver(classification_post_publish_signal, sender=Classification)
def turn_off_unsubmitted_edits(sender, classification, previously_published, newly_published, user, debug_timer: DebugTimer, **kwargs):  # pylint: disable=unused-argument
    # somewhat redundant to Classification.insert() but if publish is called without an insert, the outstanding edits flag wont be closed
    classification.flag_collection_safe.ensure_resolution(
        flag_type=classification_flag_types.classification_outstanding_edits,
        resolution='closed'
    )


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender, classification, previously_published, newly_published, user, debug_timer: DebugTimer, **kwargs):  # pylint: disable=unused-argument

    if classification.share_level_enum.index > ShareLevel.INSTITUTION.index:
        classification.flag_collection_safe.close_open_flags_of_type(
            flag_type=classification_flag_types.unshared_flag
        )

    if previously_published and previously_published.share_level in ShareLevel.DISCORDANT_LEVEL_KEYS:

        old_classification = previously_published.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        new_classification = newly_published.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)

        significant_change = ClinicalSignificance.is_significant_change(
            old_classification=old_classification,
            new_classification=new_classification
        )

        if significant_change:
            e_key = EvidenceKey.objects.get(pk=SpecialEKeys.CLINICAL_SIGNIFICANCE)
            old_label = e_key.pretty_value(old_classification)
            new_label = e_key.pretty_value(new_classification)
            # we have changed clinical significance since our last publishing
            classification.flag_collection_safe.add_flag(
                flag_type=classification_flag_types.significance_change,
                comment=f'Classification changed from {old_label} to {new_label}'
            )

            close_message = f"Classification was changed to {new_label}"

            if pending_change_flag := classification.flag_collection_safe.get_flag_of_type(classification_flag_types.classification_pending_changes):
                if pending_change_value := pending_change_flag.data.get(ClassificationFlagTypes.CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY):

                    class_url = get_url_from_view_path(
                        reverse('view_classification', kwargs={'classification_id': classification.id}),
                    )

                    pending_change_label = e_key.pretty_value(pending_change_value)
                    if new_classification != pending_change_value:
                        close_message += f", expected {pending_change_label}"

                    # is_agreed_change = new_classification != pending_change_value
                    nb = NotificationBuilder("Pending Change")
                    nb.add_markdown(
                        f"Classification <{class_url}|{classification.friendly_label}> had a pending change and has now been updated:\n"
                        f"*{old_label}* - old clinical significance\n"
                        f"*{pending_change_label}* - pending change\n"
                        f"*{new_label}* - new clinical significance")
                    nb.send()

            # note we don't care if the clin sig was was changed to agreed upon value per the flag
            # just that it was changed
            classification.flag_collection_safe.close_open_flags_of_type(
                flag_type=classification_flag_types.classification_pending_changes,
                comment=close_message
            )

    debug_timer.tick("Update share flags")
