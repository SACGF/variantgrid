from typing import Optional

from django.conf import settings
from django.dispatch import receiver
from django.urls import reverse

from classification.models.classification_import_run import (
    ClassificationImportRun,
    ClassificationImportRunStatus,
)
from library.integration_status import (
    IntegrationDetail,
    IntegrationStatus,
    integration_status_signal,
)


@receiver(signal=integration_status_signal)
def classification_import_integration_status(sender, **kwargs) -> Optional[IntegrationStatus]:
    """ Labs import on their own cadence, so no warning_age - quiet isn't broken here """
    if not settings.URLS_APP_REGISTER["classification"]:
        return None

    last_import = ClassificationImportRun.objects.order_by('-created').values_list('created', flat=True).first()
    details = [IntegrationDetail(label="Last Import", timestamp=last_import)]

    if ongoing := ClassificationImportRun.objects.filter(status=ClassificationImportRunStatus.ONGOING).count():
        details.append(IntegrationDetail(label="Ongoing Imports", extra=str(ongoing), status="warning"))

    return IntegrationStatus(
        name="Classification Imports",
        details=details,
        description="Classification records sent to us by labs",
        url=reverse('admin:classification_classificationimportrun_changelist'),
        record_count=ClassificationImportRun.objects.count()
    )
