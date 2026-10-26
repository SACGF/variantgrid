from django.db import transaction
from django.dispatch.dispatcher import receiver

from classification.enums import ShareLevel, SpecialEKeys
from classification.models.classification import (
    Classification,
    ClassificationModification,
    classification_current_state_signal,
)
from classification.models.classification_utils import ValidationMerger
from classification.models.evidence_key import EvidenceKey

INTERNAL_REVIEW_RELEVANT_DAYS = 365


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
