import django.dispatch
from django.conf import settings

from annotation.models import GeneCountType
from annotation.models.models import CachedWebResource
from annotation.tasks.cached_web_resource_tasks import ClinVarCitationsWebResourceTask

annotation_run_complete_signal = django.dispatch.Signal()

# #1670: an AnnotationRun's output is being thrown away rather than imported - a losing attempt
# unwinding, or the dispatcher scrubbing a stalled run back to CREATED. Sent with annotation_run; the
# receiver (annotation.signals.annotation_run_cleanup) reclaims the disk. Separate from
# annotation_run_complete_signal because a discarded attempt must leave the run-keyed AnnotSV directory
# alone - the attempt that won is still using it.
annotation_run_discarded_signal = django.dispatch.Signal()

clinvar_citations_post_save_handler = CachedWebResource.named_handler_factory(settings.CACHED_WEB_RESOURCE_CLINVAR_CITATIONS,
                                                                              ClinVarCitationsWebResourceTask)


def gene_counts_classification_publish_handler(sender, classification, **kwargs):  # pylint: disable=unused-argument
    GeneCountType.handle_classification_change(classification)


def gene_counts_classification_withdraw_handler(sender, classification, **kwargs):  # pylint: disable=unused-argument
    GeneCountType.handle_classification_change(classification)
