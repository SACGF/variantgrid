from functools import cached_property
from typing import Optional

import django
from attr import dataclass
from django.dispatch import Signal
from django.http import HttpResponse, JsonResponse

from library.log_utils import report_message

preview_request_signal: Signal = Signal()
"""
Receive (and return a PreviewResponse) if your module is capable of providing some tooltip previews of links
"""

@dataclass
class PreviewResponse:
    title: str
    summary: Optional[str] = None
    internal_url: Optional[str] = None
    external_url: Optional[str] = None

    def as_json(self):
        return {
            "title": self.title,
            "summary": self.summary,
            "internalUrl": self.internal_url,
            "externalUrl": self.external_url
        }


class PreviewRequest:
    """
    Similar to a SearchRequest, except everything we preview should already be formatted properly in a
    db: Database (e.g. OMIM, MONDO, PMID) and
    idx: The id within the database e.g. 453432
    
    """

    def __init__(self, db: str, idx: str):
        self.db = db
        self.idx = idx

    def preview_data(self) -> Optional[PreviewResponse]:
        results = []
        for caller, result in preview_request_signal.send_robust(sender=None, preview_request=self):
            if isinstance(result, Exception):
                report_message(f"Exception generating preview by {caller}: {result}")
            elif result is not None:
                results.append(result)
        if result_length := len(results):
            if result_length > 1:
                report_message(f"Preview data request for {self.db}:{self.idx} resulted in {result_length} responses, expected 0 or 1")
            return results[0]


def preview_view(request, db: str, idx: str):
    if preview_data := PreviewRequest(db, idx).preview_data():
        return JsonResponse(preview_data.as_json())
    else:
        return JsonResponse({"found": False})