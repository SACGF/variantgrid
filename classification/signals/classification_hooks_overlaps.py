from django.db.models.signals import post_save
from django.dispatch import receiver
from classification.models import ClassificationGrouping
from classification.services.overlaps_services import OverlapServices


@receiver(post_save, sender=ClassificationGrouping)
def classification_grouping_post_save(sender, instance: ClassificationGrouping, **kwargs):
    # a dirty grouping's summary is stale, update() saves it again once it's recalculated
    if instance.dirty:
        return
    OverlapServices.update_classification_grouping_overlap_contribution(instance)
