from typing import Optional, Union
from attr import dataclass
from django.db.models import Model
from django.dispatch import Signal
from django.http import JsonResponse
from django.utils.safestring import SafeString
from library.log_utils import report_message
from library.utils import pretty_label

preview_request_signal: Signal = Signal()
"""
Receive (and return a PreviewData) if your module is capable of providing some tooltip previews of links
"""

@dataclass
class PreviewData:
    category: str
    identifier: str
    title: Optional[str] = None
    icon: Optional[str] = None
    summary: Optional[Union[str, SafeString]] = None
    internal_url: Optional[str] = None
    external_url: Optional[str] = None

    @staticmethod
    def for_object(
            obj: Model,
            category: Optional[str] = None,
            identifier: Optional[str] = None,
            title: Optional[str] = None,
            icon: Optional[str] = None,
            summary: Optional[Union[str, SafeString]] = None,
            internal_url: Optional[str] = None,
            external_url: Optional[str] = None
        ):
        category = category or pretty_label(obj._meta.verbose_name)
        if not identifier:
            if isinstance(obj.pk, str):
                identifier = obj.pk
                if title is None:
                    title = str(obj)
            else:
                identifier = str(obj)

        if hasattr(obj, "get_absolute_url"):
            internal_url = internal_url or obj.get_absolute_url()
        return PreviewData(
            category=category,
            icon=icon,
            identifier=identifier,
            title=title,
            summary=summary,
            internal_url=internal_url,
            external_url=external_url
        )



    def as_json(self):
        summary = self.summary
        if summary and len(summary) > 200:
            summary = summary[0:200] + "..."

        return {
            "category": self.category,
            "identifier": self.identifier,
            "title": self.title,
            "icon": self.icon,
            "summary": summary,
            "internal_url": self.internal_url,
            "external_url": self.external_url
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

    def preview_data(self) -> Optional[PreviewData]:
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