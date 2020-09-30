from django.conf import settings
from django.db.models.signals import post_save
from django.dispatch.dispatcher import receiver

from variantclassification.enums import ClinicalSignificanceComparison
from variantclassification.enums.variant_classification_enums import ShareLevel, \
    SpecialEKeys
from variantclassification.models.evidence_key import EvidenceKey
from variantclassification.models.flag_types import variant_classification_flag_types
from variantclassification.models.variant_classification import \
    VariantClassification, variant_classification_post_publish_signal, \
    variant_classification_revalidate_signal


@receiver(variant_classification_revalidate_signal, sender=VariantClassification)
def revalidate(sender, variant_classification: VariantClassification, **kwargs):
    if settings.UNSHARED_FLAG_ENABLED and variant_classification.share_level_enum.index <= ShareLevel.INSTITUTION.index:
        variant_classification.flag_collection_safe.get_or_create_open_flag_of_type(
            flag_type=variant_classification_flag_types.unshared_flag,
        )


@receiver(post_save, sender=VariantClassification)
def classification_created(sender, instance, created, raw, using, update_fields, **kwargs):
    if created and settings.UNSHARED_FLAG_ENABLED:
        instance.flag_collection_safe.get_or_create_open_flag_of_type(
            flag_type=variant_classification_flag_types.unshared_flag,
        )


@receiver(variant_classification_post_publish_signal, sender=VariantClassification)
def published(sender, variant_classification, previously_published, newly_published, user, **kwargs):

    if variant_classification.share_level_enum.index > ShareLevel.INSTITUTION.index:
        variant_classification.flag_collection_safe.close_open_flags_of_type(
            flag_type=variant_classification_flag_types.unshared_flag
        )

    if previously_published and previously_published.share_level in ShareLevel.DISCORDANT_LEVEL_KEYS:

        old_classification = previously_published.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        new_classification = newly_published.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)

        significant_change = ClinicalSignificanceComparison.is_significant_change(
            old_classification=old_classification,
            new_classification=new_classification
        )

        if significant_change:
            ekey = EvidenceKey.objects.get(pk=SpecialEKeys.CLINICAL_SIGNIFICANCE)
            old_label = ekey.pretty_value(old_classification)
            new_label = ekey.pretty_value(new_classification)
            # we have changed clinical significance since our last publishing
            variant_classification.flag_collection_safe.add_flag(
                flag_type=variant_classification_flag_types.significance_change,
                comment=f'Classification changed from {old_label} to {new_label}'
            )
