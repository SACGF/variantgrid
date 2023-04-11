from typing import List

from django.dispatch import receiver

from annotation.models.models_citations import CitationIdNormalized
from library.preview_request import preview_request_signal, PreviewRequest, PreviewResponse


@receiver(signal=preview_request_signal)
def citation_preview(sender, preview_request: PreviewRequest, **kwargs):
    if preview_request.db in {"PUBMED", "PMID", "PMC", "NCBIBOOKSHELF"}:
        citation_id = CitationIdNormalized.normalize_id(preview_request.db + ":" + preview_request.idx)
        citation = citation_id.get_or_create()

        text_segments: List[str] = []

        if author_short := citation.authors_short:
            text_segments.append(author_short)
            if not citation.single_author:
                text_segments.append('et al')

        if year := citation.year:
            text_segments.append(year)

        title = citation.title
        if not title and text_segments:
            title = "Could not load title"

        text_segments.append(title)

        summary = citation.abstract or ""
        if len(summary) > 200:
            summary = summary[0:200] + "..."

        return PreviewResponse(
            title=" ".join(text_segments),
            summary=summary,
            internal_url=citation.get_absolute_url()
        )
