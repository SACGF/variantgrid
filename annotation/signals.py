from django.conf import settings
import django.dispatch

from annotation.models import GeneCountType
from annotation.models.models import CachedWebResource
from annotation.tasks.cached_web_resource_tasks import ClinGenValidityCurationsWebResourceTask

annotation_run_complete_signal = django.dispatch.Signal()

clingen_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_CLINGEN_DISEASE_VALIDITY,
                                                                    ClinGenValidityCurationsWebResourceTask)


def gene_counts_classification_publish_handler(sender, classification, **kwargs):
    GeneCountType.handle_classification_change(classification)


def gene_counts_classification_withdraw_handler(sender, classification, **kwargs):
    GeneCountType.handle_classification_change(classification)
