from django.db.models.signals import post_save
from django.dispatch import receiver

from classification.enums import ConflictType
from classification.models import ConflictLab, ClassificationGrouping, classification_grouping_clin_sig_signal, \
    classification_grouping_onc_path_signal, DiscordanceReportTriageStatus, ConflictComment
from classification.services.conflict_services import apply_conflict_lab_to_grouping, conflict_lab_for_grouping
from library.guardian_utils import admin_bot


@receiver(post_save, sender=ConflictLab)
def conflict_lab_save_update_grouping(sender, instance: ConflictLab, **kwargs):
    # Updates a ConflictGroup's will_amend values
    # so we can quickly render a warning
    apply_conflict_lab_to_grouping(instance)


@receiver(classification_grouping_clin_sig_signal, sender=ClassificationGrouping)
def clin_sig_changed(sender, instance: ClassificationGrouping, **kwargs):
    value_changed(instance, ConflictType.CLIN_SIG)

@receiver(classification_grouping_onc_path_signal, sender=ClassificationGrouping)
def onc_path_changed(sender, instance: ClassificationGrouping, **kwargs):
    value_changed(instance, ConflictType.ONCPATH)


def value_changed(instance: ClassificationGrouping, conflict_type: ConflictType):
    if conflict_lab := conflict_lab_for_grouping(instance, conflict_type):
        if conflict_lab.status == DiscordanceReportTriageStatus.REVIEWED_WILL_FIX:
            conflict_lab.status = DiscordanceReportTriageStatus.PENDING
            # note that the save will trigger an update back to the Conflict record itself
            # so this is a bit circular.
            conflict_lab.save()
            ConflictComment(
                conflict=conflict_lab.conflict,
                user=admin_bot(),
                comment="Record updated - auto changing triage status",
                meta_data={conflict_lab.lab_id: DiscordanceReportTriageStatus.PENDING}
            ).save()