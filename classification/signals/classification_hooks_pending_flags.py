from django.dispatch import receiver

from classification.models import classification_flag_types, Classification, ClassificationSummaryCalculator, \
    ClassificationModification
from flags.models import flag_comment_action, Flag, FlagComment, FlagResolution


@receiver(flag_comment_action, sender=Flag)
def check_for_pending(sender, flag_comment: FlagComment, old_resolution: FlagResolution, **kwargs):
    if flag_comment.flag.flag_type == classification_flag_types.classification_pending_changes:
        if c := Classification.objects.filter(flag_collection__id=flag_comment.flag.collection_id).first():
            print("UPDATING SUMMARY")
            c.summary = ClassificationSummaryCalculator(ClassificationModification.objects.filter(classification=c, is_lab_published=True)).cache_dict()
            c.save(update_fields=["summary"])