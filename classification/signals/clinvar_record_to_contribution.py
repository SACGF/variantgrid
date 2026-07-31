from django.dispatch.dispatcher import receiver

from annotation.models import clinvar_record_collection_refreshed, ClinVarRecordCollection
from classification.services.overlaps_services import OverlapServices


@receiver(signal=clinvar_record_collection_refreshed, sender=ClinVarRecordCollection)
def clinvar_record_collection_refreshed_to_contribution(sender, instance: ClinVarRecordCollection, **kwargs):
    if instance.expert_panel:
        OverlapServices.update_clinvar_overlap_contribution(instance)