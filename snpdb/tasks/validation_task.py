import celery

from library.log_utils import report_message
from snpdb.models.models_variant import Allele
from variantopedia.views import notify_server_status


@celery.task(track_started=True)
def validate_variant_data():
    report_message("Running automated allele validation", level="info")
    allele: Allele
    for allele in Allele.objects.all():
        allele.validate()

    notify_server_status()
