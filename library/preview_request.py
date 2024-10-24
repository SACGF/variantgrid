from dataclasses import dataclass
from datetime import datetime
from functools import cached_property
from typing import Optional, Union, Any, Callable, Type

from django.db.models import Model
from django.dispatch import Signal
from django.http import JsonResponse
from django.urls import NoReverseMatch
from django.utils.safestring import SafeString
from threadlocals.threadlocals import get_current_user

from library.log_utils import report_message
from library.utils import pretty_label
from variantgrid.perm_path import get_visible_url_names

preview_request_signal: Signal = Signal()
"""
Receive (and return a PreviewData) if your module is capable of providing some tooltip previews of links
"""


preview_extra_signal = Signal()
"""
Apps can receive this to return PreviewKeyValues for other app models, e.g. Classification can listen to this
to provide Classification counts to Alleles
"""


@dataclass
class PreviewKeyValue:
    key: Optional[str] = None
    value: Optional[Any] = None  # really we always want a value, but key is optional and want to order it first
    dedicated_row: bool = False
    icon: Optional[str] = None
    link: Optional[str] = None
    """
    Highlights this bit of data as being linked to another data set (as opposed to a KeyValue that every record
    of this type would hve. Worth highlighting)
    """

    @property
    def is_count(self):
        return isinstance(self.value, int)

    @staticmethod
    def count(preview_coordinator: 'PreviewModelMixin', amount: int, override_label: Optional[str] = None) -> Optional['PreviewKeyValue']:
        if preview_coordinator.preview_enabled():
            parts = override_label if override_label else f"{preview_coordinator.preview_category()} Count"
            return PreviewKeyValue(
                key=parts,
                value=amount,
                icon=preview_coordinator.preview_icon()
            )

    @property
    def value_str(self) -> Union[str, SafeString]:
        if isinstance(self.value, (str, SafeString)):
            return self.value
        elif isinstance(self.value, datetime):
            return f"{datetime:%Y-%m-%d}"
        else:
            return str(self.value)

    def __str__(self):
        if self.key:
            return f"{self.key}: {self.value}"
        else:
            return self.value


class PreviewProxyModel:

    def __init__(self, preview_category: str, preview_icon: str, enabled_check: Optional[Callable[[], bool]] = None):
        self._preview_category = preview_category
        self._preview_icon = preview_icon
        self._enabled_check = enabled_check

    def preview_category(self):
        return self._preview_category

    def preview_icon(self):
        return self._preview_icon

    def preview_enabled(self):
        if self._enabled_check:
            return self._enabled_check()
        else:
            return True


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
        if hasattr(cls, '_meta'):
            return pretty_label(cls._meta.verbose_name)
        else:
            return pretty_label(cls.__name__)

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
            summary_extra: Optional[list[PreviewKeyValue]] = None,
            internal_url: Optional[str] = None,
            external_url: Optional[str] = None,
            genome_builds: Optional[set['GenomeBuild']] = None,
            annotation_consortia: Optional[set['AnnotationConsortium']] = None
    ) -> 'PreviewData':
        if summary:
            if not summary_extra:
                summary_extra = []
            summary_extra.append(PreviewKeyValue(value=summary, dedicated_row=True))

        """
        Utility function to provide PreviewData for an instance with all the defaults
        provided by the mixin and the ones that PreviewData.for_object calculate.
        """
        return PreviewData.for_object(
            obj=self,
            category=category or self.preview_category(),
            identifier=identifier,
            title=title,
            icon=icon or self.preview_icon(),
            summary_extra=summary_extra,
            internal_url=internal_url,
            external_url=external_url,
            genome_builds=genome_builds,
            annotation_consortia=annotation_consortia
        )

    @property
    def preview(self) -> 'PreviewData':
        return self.preview_with()


PreviewCoordinator = Union[Type[PreviewModelMixin], PreviewProxyModel]


@dataclass(eq=True)
class PreviewData:
    category: str
    identifier: str
    title: Optional[str] = None
    icon: Optional[str] = None
    summary_extra: Optional[list[PreviewKeyValue]] = None
    internal_url: Optional[str] = None
    external_url: Optional[str] = None
    genome_builds: Optional[set['GenomeBuild']] = None
    annotation_consortia: Optional[set['AnnotationConsortium']] = None
    obj: Optional[Any] = None
    is_operation: bool = False  # indicates that the preview data is the preview of an operation to create the data
    is_error: bool = False
    # If you add any fields, be sure to modify __hash__ below...

    @staticmethod
    def for_object(
            obj: Union[Model, PreviewModelMixin],
            category: Optional[str] = None,
            identifier: Optional[str] = None,
            title: Optional[str] = None,
            icon: Optional[str] = None,
            summary_extra: Optional[list[PreviewKeyValue]] = None,
            internal_url: Optional[str] = None,
            external_url: Optional[str] = None,
            genome_builds: Optional[set['GenomeBuild']] = None,
            annotation_consortia: Optional[set['AnnotationConsortium']] = None,
            is_operation: bool = False,
            is_error: bool = False):

        if category is None:
            if hasattr(obj, "_meta"):
                category = obj._meta.verbose_name
            else:
                category = pretty_label(obj.__class__.__name__)

        if identifier is None:
            if hasattr(obj, "pk") and isinstance(obj.pk, str):
                identifier = obj.pk
                if title is None:
                    test_title = str(obj)
                    if test_title != identifier:
                        title = test_title
            else:
                identifier = str(obj)

        if internal_url is None and hasattr(obj, "get_absolute_url"):
            try:
                internal_url = obj.get_absolute_url()
            except NoReverseMatch:
                internal_url = "javascript:alert('Could not load a view for this type of result')"

        if genome_builds is None:
            if hasattr(obj, "genome_build") and (genome_build := obj.genome_build):
                genome_builds = {genome_build}
            elif hasattr(obj, "genome_builds"):
                genome_builds = obj.genome_builds

        if annotation_consortia is None:
            if hasattr(obj, "annotation_consortium"):
                annotation_consortia = {obj.annotation_consortium}
            elif hasattr(obj, "annotation_consortia"):
                annotation_consortia = obj.annotation_consortia

        # make sure we're dealing with AnnotationConsortium and not pure string (easy to do when dealing with TextChoices)
        if annotation_consortia:
            def convert_consortium(obj_c: Any):
                if isinstance(obj_c, str):
                    from genes.models_enums import AnnotationConsortium
                    return AnnotationConsortium(obj_c)
                else:
                    return obj_c

            annotation_consortia = {convert_consortium(obj) for obj in annotation_consortia}

        return PreviewData(
            obj=obj,
            category=category,
            icon=icon,
            identifier=identifier,
            title=title,
            summary_extra=summary_extra,
            internal_url=internal_url,
            external_url=external_url,
            genome_builds=genome_builds,
            annotation_consortia=annotation_consortia,
            is_operation=is_operation,
            is_error=is_error,
        )

    @cached_property
    def summary_all(self):
        return (self.summary_extra or []) + self._summary_external()

    def _summary_external(self):
        external_extra = []
        if obj := self.obj:
            # TODO getting current user isn't great, maybe PreviewData should have a reference to the user
            user = get_current_user()
            for caller, response in preview_extra_signal.send_robust(sender=obj.__class__, user=user, obj=obj):
                if response:
                    if isinstance(response, Exception):
                        report_message("Error during preview", 'error',
                                       extra_data={"target": str(response), "caller": str(caller)})
                    else:
                        external_extra.extend(response)
        return external_extra

    @property
    def html_element(self) -> str:
        if self.is_error:
            return "span"
        return "a"

    @property
    def summary(self):
        if extra := self.summary_all:
            return ", ".join(str(kv) for kv in extra)
        return None

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

    def __hash__(self):
        genome_builds = []
        if self.genome_builds:
            genome_builds = tuple(sorted(self.genome_builds))
        annotation_consortia = []
        if self.annotation_consortia:
            annotation_consortia = tuple(sorted(self.annotation_consortia))

        fields = (
            self.category,
            self.identifier,
            self.title,
            self.icon,
            self.summary_extra,
            self.internal_url,
            self.external_url,
            genome_builds,
            annotation_consortia,
            self.obj,
            self.is_operation,
        )
        return hash((f for f in fields if f is not None))

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


def preview_view(request, db: str, idx: str) -> JsonResponse:
    if preview_data := PreviewRequest(db, idx).preview_data():
        return JsonResponse(preview_data.as_json())
    else:
        return JsonResponse({"found": False})
