from django.contrib.auth.models import User
from django.db.models.signals import pre_delete
from django.dispatch import receiver

from classification.enums import ShareLevel, SubmissionSource
from classification.models import classification_post_publish_signal, Classification, ClassificationModification, \
    classification_withdraw_signal, allele_info_changed_signal, ImportedAlleleInfo, ClassificationImportRun
from classification.models.classification_grouping import ClassificationGrouping, ClassificationGroupingEntry
from classification.models.classification_import_run import classification_imports_complete_signal
from library.utils import DebugTimer

###
# Capture events that could require a ClassificationGrouping to update.
# If the ClassificationGrouping for a classificaiton might change (or be created for the first time)
# call assign_grouping_for_classification, if it already exists and we just need to update it
# call dirty_up
##

def _instant_undirty_check():
    if ClassificationImportRun.ongoing_imports():
        pass
    else:
        ClassificationGrouping.update_all_dirty()


@receiver(classification_imports_complete_signal, sender=ClassificationImportRun)
def classification_imports_complete(sender, **kwargs):
    ClassificationGrouping.update_all_dirty()


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender,
              classification: Classification,
              previously_published: ClassificationModification,
              newly_published: ClassificationModification,
              previous_share_level: ShareLevel,
              user: User,
              debug_timer: DebugTimer,
              **kwargs):  # pylint: disable=unused-argument
    # CLASSIFICATION PUBLISHED
    # we call assign_grouping_for_classification as we might have a new allele origin,
    # this might be the first time we're publishing... so dirty_up
    ClassificationGrouping.assign_grouping_for_classification(classification)
    _instant_undirty_check()


@receiver(pre_delete, sender=Classification)
def deleting_classification(sender, instance: Classification, **kwargs):  # pylint: disable=unused-argument
    # CLASSIFICATION DELETED
    # when a classification is deleted, it will delete the corresponding ClassificationGroupingEntry
    # so we need to mark the Allele Origin grouping (which dirty up does)
    if entry := ClassificationGroupingEntry.objects.filter(classification=instance).first():
        entry.dirty_up()


@receiver(classification_withdraw_signal, sender=Classification)
def withdraw_changed(sender, classification: Classification, **kwargs):  # pylint: disable=unused-argument
    # WITHDRAW
    if entry := ClassificationGroupingEntry.objects.filter(classification=classification).first():
        entry.dirty_up()
        _instant_undirty_check()


@receiver(allele_info_changed_signal, sender=ImportedAlleleInfo)
def allele_info_changed_update_classifications(sender, allele_info: ImportedAlleleInfo, **kwargs):
    # CLASSIFICATION MATCHED TO AN ALLELE
    for classification in allele_info.classification_set.all():
        ClassificationGrouping.assign_grouping_for_classification(classification)
    _instant_undirty_check()
