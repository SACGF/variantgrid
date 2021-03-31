from django.conf import settings
from django.db.models.signals import post_save
from django.dispatch.dispatcher import receiver

from classification.enums.classification_enums import ShareLevel, SpecialEKeys, ClinicalSignificance
from classification.models.evidence_key import EvidenceKey
from classification.models.flag_types import classification_flag_types
from classification.models.classification import \
    Classification, classification_post_publish_signal, \
    classification_revalidate_signal


@receiver(classification_revalidate_signal, sender=Classification)
def revalidate(sender, classification: Classification, **kwargs):
    if settings.UNSHARED_FLAG_ENABLED and classification.share_level_enum.index <= ShareLevel.INSTITUTION.index:
        classification.flag_collection_safe.get_or_create_open_flag_of_type(
            flag_type=classification_flag_types.unshared_flag,
        )


@receiver(post_save, sender=Classification)
def classification_created(sender, instance, created, raw, using, update_fields, **kwargs):
    if created and settings.UNSHARED_FLAG_ENABLED:
        instance.flag_collection_safe.get_or_create_open_flag_of_type(
            flag_type=classification_flag_types.unshared_flag,
        )


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender, classification, previously_published, newly_published, user, **kwargs):

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
            ekey = EvidenceKey.objects.get(pk=SpecialEKeys.CLINICAL_SIGNIFICANCE)
            old_label = ekey.pretty_value(old_classification)
            new_label = ekey.pretty_value(new_classification)
            # we have changed clinical significance since our last publishing
            classification.flag_collection_safe.add_flag(
                flag_type=classification_flag_types.significance_change,
                comment=f'Classification changed from {old_label} to {new_label}'
            )
