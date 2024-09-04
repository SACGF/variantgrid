from django.contrib.auth.models import User
from django.db.models.signals import pre_delete
from django.dispatch import receiver

from classification.enums import ShareLevel
from classification.models import classification_post_publish_signal, Classification, ClassificationModification, \
    classification_withdraw_signal
from classification.models.classification_grouping import ClassificationGrouping, ClassificationGroupingEntry
from library.utils import DebugTimer


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender,
              classification: Classification,
              previously_published: ClassificationModification,
              newly_published: ClassificationModification,
              previous_share_level: ShareLevel,
              user: User,
              debug_timer: DebugTimer,
              **kwargs):  # pylint: disable=unused-argument
    ClassificationGrouping.assign_grouping_for_classification(classification)


@receiver(pre_delete, sender=Classification)
def deleted_variant(sender, instance: Classification, **kwargs):  # pylint: disable=unused-argument
    # when a classification is deleted, it will delete the corresponding ClassificationGroupingEntry
    # so we need to mark the overall grouping
    if entry := ClassificationGroupingEntry.objects.filter(classification=instance).first():
        entry.dirty_up()


@receiver(classification_withdraw_signal, sender=Classification)
def withdraw_changed(sender, classification: Classification, **kwargs):  # pylint: disable=unused-argument
    if entry := ClassificationGroupingEntry.objects.filter(classification=classification).first():
        entry.dirty_up()
