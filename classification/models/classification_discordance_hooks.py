from typing import Set, Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.db import transaction
from django.db.models.signals import post_delete
from django.dispatch.dispatcher import receiver

from classification.enums import SpecialEKeys, ShareLevel
from classification.models.classification import Classification, \
    classification_current_state_signal, \
    classification_post_publish_signal, \
    classification_variant_set_signal, ClassificationModification, \
    classification_withdraw_signal
from classification.models.classification_utils import ValidationMerger
from classification.models.clinical_context_models import ClinicalContext, \
    clinical_context_signal, ClinicalContextRecalcTrigger
from classification.models.discordance_models import DiscordanceReport
from classification.models.evidence_key import EvidenceKey, EvidenceKeyMap
from classification.models.flag_types import classification_flag_types
from library.utils import DebugTimer

INTERNAL_REVIEW_RELEVANT_DAYS = 365


@receiver(classification_variant_set_signal, sender=Classification)
def variant_set(sender, **kwargs):  # pylint: disable=unused-argument
    record: Classification = kwargs['classification']
    variant = kwargs['variant']

    old_clinical_context: Optional[ClinicalContext] = None
    new_clinical_context: Optional[ClinicalContext] = None

    if not variant and record.clinical_context:
        old_clinical_context = record.clinical_context
        record.clinical_context = None
        record.save()

    elif variant and ((not record.clinical_context) or
                      (record.clinical_context and record.clinical_context.allele != variant.allele)):
        if dg := ClinicalContext.default_group_for(variant):
            old_clinical_context = record.clinical_context
            new_clinical_context = dg
            record.clinical_context = dg
            record.save()

    # make sure you do your calculations after
    if old_clinical_context:
        old_clinical_context.recalc_and_save(cause=f'Classification {record.friendly_label} unmatched', cause_code=ClinicalContextRecalcTrigger.VARIANT_SET)
    if new_clinical_context:
        new_clinical_context.recalc_and_save(cause=f'Classification {record.friendly_label} submitted', cause_code=ClinicalContextRecalcTrigger.VARIANT_SET)


@receiver(post_delete, sender=Classification)
def deleted_variant(sender, instance: Classification, **kwargs):  # pylint: disable=unused-argument
    classification = instance
    if classification.clinical_context:
        cause = f'Classification {instance.friendly_label} deleted'
        classification.clinical_context.recalc_and_save(cause=cause, cause_code=ClinicalContextRecalcTrigger.WITHDRAW_DELETE)


@receiver(classification_withdraw_signal, sender=Classification)
def withdraw_changed(sender, classification: Classification, **kwargs):  # pylint: disable=unused-argument
    if classification.clinical_context:
        cause = ''
        if classification.withdrawn:
            cause = f'Classification {classification.friendly_label} withdrawn'
        else:
            cause = f'Classification {classification.friendly_label} un-withdrawn'

        classification.clinical_context.recalc_and_save(cause=cause, cause_code=ClinicalContextRecalcTrigger.WITHDRAW_DELETE)


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender,
              classification: Classification,
              previously_published: ClassificationModification,
              newly_published: ClassificationModification,
              previous_share_level: ShareLevel,
              user: User,
              debug_timer: DebugTimer,
              **kwargs):  # pylint: disable=unused-argument
    """
    Only care about publicly shared records
    """
    new_share_level = newly_published.share_level_enum
    if not new_share_level.is_discordant_level:
        return

    change_share_level = new_share_level != previous_share_level
    first_publish = not previous_share_level.is_discordant_level
    cs = EvidenceKey.objects.get(pk=SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(
        classification.evidence.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)) or 'Unclassified'

    diff_keys: Set[str] = set()
    if not first_publish:
        if previously_published:
            keys = set(newly_published.evidence) | set(previously_published.evidence)
            for key in keys:
                if newly_published.get(key) != previously_published.get(key):
                    diff_keys.add(key)

    cs_changed = SpecialEKeys.CLINICAL_SIGNIFICANCE in diff_keys
    at_least_minor_diffs = previously_published != newly_published

    if first_publish or change_share_level or diff_keys or at_least_minor_diffs:
        diff_str = None
        if diff_keys:
            ekeys = EvidenceKeyMap.instance()
            diff_labels = [ekeys.get(key).pretty_label for key in diff_keys]
            diff_labels.sort()
            diff_str = ', '.join(diff_labels)

        comment_action = 'Submitted' if first_publish or change_share_level else 'Resubmitted'
        comment_as = f'as {cs}' if first_publish or cs_changed else None
        comment_changed = f'with changes to {diff_str}' if diff_str and not first_publish else None
        comment_to = f'to {new_share_level.context_label(classification)}' if change_share_level else None

        comment = ' '.join([part for part in [comment_action, comment_as, comment_changed, comment_to] if part])

        # worth reporting the flag
        # don't want to report a flag if we're not making any changes at all
        classification.flag_collection_safe.add_flag(
            flag_type=classification_flag_types.submitted_flag,
            user=user,
            permission_check=False,
            comment=comment
        )

    if classification.variant and not classification.clinical_context:
        classification.clinical_context = ClinicalContext.default_group_for(classification.variant)
        classification.save(update_fields=['clinical_context'])

    if classification.clinical_context:
        cause = ''
        if first_publish:
            cause = f'Classification {classification.friendly_label} submitted'
        else:
            cause = f'Classification {classification.friendly_label} re-submitted as {cs}'

        classification.clinical_context.recalc_and_save(cause=cause, cause_code=ClinicalContextRecalcTrigger.SUBMISSION)

    debug_timer.tick("Update Clinical Grouping")


@receiver(clinical_context_signal, sender=ClinicalContext)
def clinical_context_update(sender, clinical_context: ClinicalContext, status: str, is_significance_change: bool, cause: str, **kwargs):  # pylint: disable=unused-argument
    if settings.DISCORDANCE_ENABLED:

        discordant_classifications = set()
        latest = DiscordanceReport.update_latest(clinical_context, cause=cause, update_flags=True)

    else:
        clinical_context.flag_collection_safe.close_open_flags_of_type(
            flag_type=classification_flag_types.clinical_context_discordance,
            comment='Discordance functionality has been disabled'
        )
        for vc in Classification.objects.filter(clinical_context=clinical_context):
            vc.flag_collection_safe.close_open_flags_of_type(
                flag_type=classification_flag_types.discordant,
                comment='Discordance functionality has been disabled'
            )


@transaction.atomic
@receiver(classification_current_state_signal, sender=Classification)
def discordance_current_state(sender, **kwargs) -> ValidationMerger:  # pylint: disable=unused-argument
    """
    Reports if the classification is in discordance or review
    """
    record = kwargs.get('record')  # type: Classification
    user = kwargs.get('user')  # type: User

    messages = ValidationMerger()
    my_classification = record.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)

    if record.clinical_context and my_classification:
        cs_key = EvidenceKey.objects.get(pk=SpecialEKeys.CLINICAL_SIGNIFICANCE)
        qs = ClassificationModification.latest_for_user(
            published=True,
            user=user,
            clinical_context=record.clinical_context_id,
            exclude_withdrawn=False)
        qs = qs.exclude(classification=record)

        different_vcs = []

        for vcm in qs:
            other_cs = vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            if other_cs:
                is_same = other_cs == my_classification
                if not is_same:
                    other_record = vcm.classification
                    other_title = other_record.lab.name + ' / ' + other_record.lab_record_id

                    if record.share_level not in ShareLevel.DISCORDANT_LEVEL_KEYS or \
                        vcm.share_level not in ShareLevel.DISCORDANT_LEVEL_KEYS:

                        different_vcs.append(vcm.classification)

        if different_vcs:
            message = 'Different clinical significance compared to\n'
            for other_record in different_vcs:
                other_title = other_record.lab.name + ' / ' + other_record.lab_record_id
                other_cs = other_record.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
                other_classification = cs_key.pretty_value(other_cs)
                message = message + f'{other_title} - {other_classification}\n'

            messages.add_message(
                key=SpecialEKeys.CLINICAL_SIGNIFICANCE,
                severity='info',
                code='private-discordance',
                message=message,
                link=record.variant.get_absolute_url()
            )

    return messages
