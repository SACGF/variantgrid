from django.contrib.auth.models import User
from django.dispatch.dispatcher import receiver

from classification.models import ClinicalContext
from classification.models.classification_import_run import ClassificationImportRun, \
    classification_imports_complete_signal
from classification.models.clinical_context_models import ClinicalContextRecalcTrigger
from flags.models.models import FlagCollection, \
    flag_collection_extra_info_signal, FlagInfos
from classification.signals import send_prepared_discordance_notifications


@receiver(flag_collection_extra_info_signal, sender=FlagCollection)
def get_extra_info(flag_infos: FlagInfos, user: User, **kwargs):  # pylint: disable=unused-argument
    ccs = ClinicalContext.objects.filter(flag_collection__in=flag_infos.ids).select_related('allele')
    for cc in ccs:
        flag_infos.set_extra_info(cc.flag_collection_id, {
            'label': f'{cc.allele} - "{cc.name}"',
            'cc_id': cc.id
        }, source_object=cc)


@receiver(classification_imports_complete_signal, sender=ClassificationImportRun)
def import_complete(**kwargs):
    # this is called when there are no ongoing imports, find all the delayed clinical contexts, and calculate them
    for cc in ClinicalContext.objects.filter(pending_cause__isnull=False):
        # cause should automatically be loaded from pending cause anyway
        cc.recalc_and_save(cause=cc.pending_cause, cause_code=ClinicalContextRecalcTrigger.DELAYED)

    send_prepared_discordance_notifications()
