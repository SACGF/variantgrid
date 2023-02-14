from django.conf import settings
from django.db.models.signals import post_save
from django.dispatch.dispatcher import receiver

from classification.enums.classification_enums import ShareLevel
from classification.models.classification import \
    Classification, classification_revalidate_signal
from classification.models.flag_types import classification_flag_types


@receiver(classification_revalidate_signal, sender=Classification)
def revalidate(sender, classification: Classification, **kwargs):  # pylint: disable=unused-argument
    if settings.UNSHARED_FLAG_ENABLED and classification.share_level_enum.index <= ShareLevel.INSTITUTION.index:
        classification.flag_collection_safe.get_or_create_open_flag_of_type(
            flag_type=classification_flag_types.unshared_flag,
        )


@receiver(post_save, sender=Classification)
def classification_created(sender, instance, created, raw, using, update_fields, **kwargs):  # pylint: disable=unused-argument
    if created and settings.UNSHARED_FLAG_ENABLED:
        instance.flag_collection_safe.get_or_create_open_flag_of_type(
            flag_type=classification_flag_types.unshared_flag,
        )
