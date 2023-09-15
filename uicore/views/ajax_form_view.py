from dataclasses import dataclass
from enum import Enum
from typing import TypeVar, Dict, Generic, Callable, Optional, Tuple

from django.db.models import Model
from django.http import HttpResponse, HttpRequest
from django.template import loader
from django.utils.safestring import SafeString
from django.views import View

"""
This code is designed to handle showing data about an object in a myriad of ways, but with minimal code duplication.
Supports:
Inline with popup modal that replaces the inline on submit and
Card (that is then ajax replaces itself on edit & submit)

Note the one thing it currently can't do is embed an editable form (but can embed a card view that turns into an editable form through ajax)

relies on global.js having support for:
embedded-wrapper (for replacing context of card with editable version, with fallback to modal contents)
data-replace (for replacing contents of inline if it opened up a form)
"""


T = TypeVar("T", bound=Model)


class AjaxFormMode(str, Enum):
    INLINE = "inline"
    CARD = "card"
    EMBEDDED_CARD = "embedded-card"  # The first time you embed a card, do this, so it can be replaced
    MODAL = "modal"  # almost identical to INLINE but on save will trigger the modal to disappear


@dataclass
class LazyRender(Generic[T]):
    """
    A delayed rendering, that may be delayed, so a template can include it
    or we might call render() on it right after.
    It helps in reducing the duplication of code, your AJAX
    """

    core_object: T
    core_object_name: str
    template_name: str
    """
    Has a successful save on the form just occurred
    """

    static_context: Optional[Dict] = None
    """
    A dictionary that can be passed into the constructor that will be part of the context provided to the template
    """

    dynamic_context: Optional[Callable[[HttpRequest], Dict]] = None
    """
    A callable that takes a request and produces a dictionary for context. Is combined with static context.
    """

    def _context(self, request):
        """
        Create the context by combining the static and dynamic context
        """
        use_context = (self.static_context or {}).copy()
        if dynamic_context := self.dynamic_context:
            use_context.update(dynamic_context(request))
        if self.core_object_name not in use_context:
            use_context[self.core_object_name] = self.core_object
        return use_context

    def _wrap(self, mode: AjaxFormMode, saved: bool = False) -> Tuple[str, AjaxFormMode, str]:
        """
        Produces:
        starting HTML
        the mode to send to the template
        ending HTML
        based on the mode provided and the saved status
        """
        if mode == "inline":
            return (
                f'<div class="d-inline-block" id="{self.core_object_name}-{ self.core_object.pk }">',
                AjaxFormMode.INLINE,
                '</div>'
            )
        elif mode == AjaxFormMode.MODAL and saved:
            return (
                f'<div class="d-none auto-close-modal"><div class="d-inline-block" data-replace="#{self.core_object_name}-{self.core_object.pk}">',
                AjaxFormMode.INLINE,
                '</div></div>'
            )
        elif mode == AjaxFormMode.EMBEDDED_CARD:
            return (
                f'<div class="card embed-wrapper">',
                AjaxFormMode.CARD,
                '</div>'
            )
        else:
            return (
                '<div class="card">',
                mode,  # could be Modal or Card
                '</div>'
            )

    def _render_to_string(self, request, context, mode: AjaxFormMode, saved: bool = False) -> str:
        parts = self._wrap(mode, saved=saved)
        context["mode"] = parts[1].value
        template_output = loader.render_to_string(template_name=self.template_name, context=context, request=request)
        return parts[0] + template_output + parts[2]

    def render(self, request, saved: bool = False) -> HttpResponse:
        mode = AjaxFormMode(request.GET.get("mode", "card"))
        use_context = self._context(request)
        content_str = self._render_to_string(request, use_context, mode=mode, saved=saved)
        return HttpResponse(bytes(content_str, "utf-8"))

    def embed(self, request, mode: AjaxFormMode = AjaxFormMode.EMBEDDED_CARD, **kwargs) -> SafeString:
        """
        Returns a safe string you can inject directly into a template
        """
        use_context = self._context(request)
        use_context.update(kwargs)
        return SafeString(self._render_to_string(request, use_context, mode=mode))


class AjaxFormView(View, Generic[T]):
    """
    If you want to be able to edit a form within a page - and use AJAX to submit and refresh the form, extend this.
    """

    @classmethod
    def lazy_render(cls, obj: T, context: Optional[Dict] = None) -> LazyRender:
        """
        This method can be called to give a LazyRender that can be embedded as needed
        or it can be called by this view on get/post followed by .render()
        """
        pass
