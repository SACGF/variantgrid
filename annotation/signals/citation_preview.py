from typing import List

from django.dispatch import receiver

from annotation.models import CitationFetchRequest
from annotation.models.models_citations import CitationIdNormalized
from library.preview_request import preview_request_signal, PreviewRequest, PreviewData


@receiver(signal=preview_request_signal)
def citation_preview(sender, preview_request: PreviewRequest, **kwargs):
    if preview_request.db in {"PUBMED", "PMID", "PMC", "NCBIBOOKSHELF"}:
        citation_id = CitationIdNormalized.normalize_id(preview_request.db + ":" + preview_request.idx)
        response = CitationFetchRequest.fetch_all_now([citation_id])
        return response.for_requested(citation_id).preview
