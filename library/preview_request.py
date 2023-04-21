from abc import ABC, abstractmethod
from typing import Optional, Union, Any
from attr import dataclass
from django.db.models import Model
from django.dispatch import Signal
from django.http import JsonResponse
from django.urls import NoReverseMatch
from django.utils.safestring import SafeString
from library.log_utils import report_message
from library.utils import pretty_label
from variantgrid.perm_path import get_visible_url_names

preview_request_signal: Signal = Signal()
"""
Receive (and return a PreviewData) if your module is capable of providing some tooltip previews of links
"""


class PreviewModelMixin:
    """
    Having models implement this class makes it easier for search results to give summaries, PreviewData
    can be used in other contexts like providing information about IDs found in text
    """

    @classmethod
    def preview_category(cls) -> str:
        """
        Provide an English language name for this data
        """
        return pretty_label(cls._meta.verbose_name)

    @classmethod
    def preview_if_url_visible(cls) -> Optional[str]:
        """
        If this app is included, but this model is sometimes not considered active, what URL must be visible
        for this app to surface
        """
        return None

    @classmethod
    def preview_enabled(cls) -> bool:
        """
        If this returns false, pretend like the model doesn't exist, don't search for it
        """
        if url_name := cls.preview_if_url_visible():
            return bool(get_visible_url_names().get(url_name))
        return True

    @classmethod
    def preview_icon(cls) -> str:
        """
        A set of classes (typically font awesome) to put in a <i class_name=""></i> to give a visual clue as to the data
        """
        return "fa-solid fa-circle"

    def preview_with(
            self,
            category: Optional[str] = None,
            identifier: Optional[str] = None,
            title: Optional[str] = None,
            icon: Optional[str] = None,
            summary: Optional[Union[str, SafeString]] = None,
            internal_url: Optional[str] = None,
            external_url: Optional[str] = None,
            genome_build: Optional['GenomeBuild'] = None,
            annotation_consortium: Optional['AnnotationConsortium'] = None
    ) -> 'PreviewData':

        """
        Utility function to provide PreviewData for an instance with all the defaults
        provided by the mixin and the ones that PreviewData.for_object calculate.
        """
        # TODO maybe get rid of this method and have PreviewData.for_object check all the defaults here
        return PreviewData.for_object(
            self,
            category=category or self.preview_category(),
            identifier=identifier,
            title=title,
            icon=icon or self.preview_icon(),
            summary=summary,
            internal_url=internal_url,
            external_url=external_url,
            genome_build=genome_build,
            annotation_consortium=annotation_consortium
        )

    @property
    def preview(self) -> 'PreviewData':
        return self.preview_with()


@dataclass
class PreviewData:
    category: str
    identifier: str
    title: Optional[str] = None
    icon: Optional[str] = None
    summary: Optional[Union[str, SafeString]] = None
    internal_url: Optional[str] = None
    external_url: Optional[str] = None
    genome_build: Optional['GenomeBuild'] = None
    annotation_consortium: Optional['AnnotationConsortium'] = None
    obj: Optional[Any] = None

    @staticmethod
    def for_object(
            obj: Model,
            category: Optional[str] = None,
            identifier: Optional[str] = None,
            title: Optional[str] = None,
            icon: Optional[str] = None,
            summary: Optional[Union[str, SafeString]] = None,
            internal_url: Optional[str] = None,
            external_url: Optional[str] = None,
            genome_build: Optional['GenomeBuild'] = None,
            annotation_consortium: Optional['AnnotationConsortium'] = None):

        if category is None:
            if hasattr(obj, "_meta"):
                category = obj._meta.verbose_name
            else:
                category = pretty_label(obj.__class__.__name__)

        if identifier is None:
            if hasattr(obj, "pk") and isinstance(obj.pk, str):
                identifier = obj.pk
                if title is None:
                    title = str(obj)
            else:
                identifier = str(obj)

        if internal_url is None and hasattr(obj, "get_absolute_url"):
            try:
                internal_url = obj.get_absolute_url()
            except NoReverseMatch:
                internal_url = "javascript:alert('Could not load a view for this type of result')"

        if genome_build is None and hasattr(obj, "genome_build"):
            genome_build = obj.genome_build

        if annotation_consortium is None and hasattr(obj, "annotation_consortium"):
            annotation_consortium = obj.annotation_consortium

        return PreviewData(
            obj=obj,
            category=category,
            icon=icon,
            identifier=identifier,
            title=title,
            summary=summary,
            internal_url=internal_url,
            external_url=external_url,
            genome_build=genome_build,
            annotation_consortium=annotation_consortium
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