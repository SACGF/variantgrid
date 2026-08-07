from django.dispatch import receiver

from library.preview_request import PreviewRequest, preview_request_signal
from ontology.models import OntologyTerm


@receiver(signal=preview_request_signal)
def preview_ontology(sender, preview_request: PreviewRequest, **kwargs):
    if preview_request.db in {"MONDO", "OMIM", "HPO"}:
        return OntologyTerm.get_or_stub(preview_request.db + ":" + preview_request.idx).preview
