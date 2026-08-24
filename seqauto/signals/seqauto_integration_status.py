from django.conf import settings
from django.dispatch import receiver
from django.urls import reverse

from library.integration_status import (
    IntegrationDetail,
    IntegrationStatus,
    integration_status_signal,
)
from seqauto.models import QCExecSummary, SampleSheet, SequencingRun


def _latest_created(model):
    return model.objects.order_by('-created').values_list('created', flat=True).first()


@receiver(signal=integration_status_signal)
def seqauto_integration_status(sender, **kwargs) -> list[IntegrationStatus]:
    """ The sequencing pipeline posts into /seqauto/api/ - these are the records it leaves behind.
        Sequencing happens on the lab's schedule, so no warning_age - quiet isn't broken here """
    if not settings.SEQAUTO_ENABLED:
        return []

    return [IntegrationStatus(
        name="Sequencing Pipeline",
        details=[
            IntegrationDetail(label="Last Sequencing Run", timestamp=_latest_created(SequencingRun)),
            IntegrationDetail(label="Last Sample Sheet", timestamp=_latest_created(SampleSheet)),
            IntegrationDetail(label="Last QC", timestamp=_latest_created(QCExecSummary)),
        ],
        description="Sequencing runs, sample sheets and QC sent up by SeqAuto",
        url=reverse('sequencing_runs'),
        record_count=SequencingRun.objects.count()
    )]
