from django.dispatch import receiver

from classification.models.clinvar_export_models import (
    ClinVarExportBatch,
    ClinVarExportBatchStatus,
)
from classification.models.clinvar_export_sync import clinvar_export_sync
from library.integration_status import (
    IntegrationDetail,
    IntegrationDirection,
    IntegrationStatus,
    integration_status_signal,
)
from snpdb.models import ClinVarKey


def _latest_batch(clinvar_key: ClinVarKey, status: str) -> ClinVarExportBatch:
    return ClinVarExportBatch.objects.filter(clinvar_key=clinvar_key, status=status).order_by('-created').first()


@receiver(signal=integration_status_signal)
def clinvar_export_integration_status(sender, **kwargs) -> list[IntegrationStatus]:
    if not clinvar_export_sync.is_enabled:
        return []

    statuses: list[IntegrationStatus] = []
    for clinvar_key in ClinVarKey.objects.order_by('pk'):
        last_submitted = _latest_batch(clinvar_key, ClinVarExportBatchStatus.SUBMITTED)
        last_rejected = _latest_batch(clinvar_key, ClinVarExportBatchStatus.REJECTED)
        if not (last_submitted or last_rejected):
            continue

        details = [IntegrationDetail(
            label="Last Submitted",
            timestamp=last_submitted.created if last_submitted else None
        )]
        if last_rejected:
            newer_than_submitted = not last_submitted or last_rejected.created > last_submitted.created
            details.append(IntegrationDetail(
                label="Last Rejected",
                timestamp=last_rejected.created,
                status="danger" if newer_than_submitted else "warning",
                extra=last_rejected.clinvar_batch_id
            ))

        latest_batch = last_rejected if (last_rejected and (not last_submitted or last_rejected.created > last_submitted.created)) else last_submitted
        statuses.append(IntegrationStatus(
            name=f"ClinVar Export ({clinvar_key.name or clinvar_key.pk})",
            details=details,
            direction=IntegrationDirection.OUTBOUND,
            url=latest_batch.get_absolute_url(),
            record_count=ClinVarExportBatch.objects.filter(clinvar_key=clinvar_key).count()
        ))
    return statuses
