from django.db.models.signals import post_save
from django.dispatch import receiver

from classification.models import ConflictLab
from classification.services.conflict_services import apply_conflict_lab_to_grouping


@receiver(post_save, sender=ConflictLab)
def conflict_lab_save_update_grouping(sender, instance: ConflictLab, **kwargs):
    # Updates a ConflictGroup's will_amend values
    # so we can quickly render a warning
    apply_conflict_lab_to_grouping(instance)