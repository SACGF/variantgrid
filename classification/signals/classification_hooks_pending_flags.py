from django.dispatch import receiver

from classification.models import (
    Classification,
    ClassificationSummaryCalculator,
    classification_flag_types,
)
from flags.models import Flag, FlagComment, FlagResolution, flag_comment_action


@receiver(flag_comment_action, sender=Flag)
def check_for_pending(sender, flag_comment: FlagComment, old_resolution: FlagResolution, **kwargs):
    if flag_comment.flag.flag_type == classification_flag_types.classification_pending_changes:
        if c := Classification.objects.filter(flag_collection__id=flag_comment.flag.collection_id).first():
            # print("UPDATING SUMMARY")
            c.summary = ClassificationSummaryCalculator(c.last_published_version).cache_dict()
            c.save(update_fields=["summary"])
