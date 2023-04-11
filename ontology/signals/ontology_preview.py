from django.dispatch import receiver

from library.preview_request import preview_request_signal, PreviewRequest, PreviewResponse
from ontology.models import OntologyTerm


@receiver(signal=preview_request_signal)
def preview_ontology(sender, preview_request: PreviewRequest, **kwargs):
    if preview_request.db in {"MONDO", "OMIM", "HPO"}:
        term = OntologyTerm.get_or_stub(preview_request.db + ":" + preview_request.idx)
        if term.is_stub:
            name = "Term not found"
        else:
            name = term.name

        # can introduce warning to preview response, but worried this will be overkill
        # if OMIMs turned up outside of condition field
        # warning = None
        # if not term.is_valid_for_condition:
        #     warning = "Not a valid term for condition"

        return PreviewResponse(
            title=name,
            summary=term.definition,
            internal_url=term.get_absolute_url()
        )