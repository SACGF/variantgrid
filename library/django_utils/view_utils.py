"""
View helpers: is_ajax, render_ajax_view (a template rendered bare for an AJAX load, or wrapped in
snpdb/embedded_ajax.html when opened as a page), view_to_string (a resolved view's dotted name) and
login_not_required_include (a third-party URLconf that LoginRequiredMiddleware must let through).
"""
import inspect
from importlib import import_module
from typing import Optional

from django.contrib.auth.decorators import login_not_required
from django.http import HttpRequest, HttpResponse
from django.shortcuts import render
from django.template.loader import render_to_string
from django.urls import URLPattern, URLResolver, include


def is_ajax(request: HttpRequest):
    return request.META.get('HTTP_X_REQUESTED_WITH') == 'XMLHttpRequest'


def render_ajax_view(
        request: HttpRequest,
        template_name: str,
        context: dict,
        menubar: Optional[str] = None) -> HttpResponse:
    if not context:
        context = {}
    if is_ajax(request):
        context['render_mode'] = 'ajax'
        return render(request, template_name, context)
    else:
        context['render_mode'] = 'embedded'
        text = render_to_string(template_name, context, request=request)
        return render(request, "snpdb/embedded_ajax.html", {"embedded_content": text, "menubar": menubar})


def view_to_string(view) -> str:
    """ This isn't fast so don't use this for anything other than displaying info in
        error pathways or something """
    v=inspect.unwrap(view)

    vc=getattr(v,'view_class',None)
    if vc is not None:
        m=vc.__module__
        n=vc.__name__
        return f'{m}.{n}'

    if inspect.ismethod(v):
        c=v.__self__.__class__
        return f'{c.__module__}.{c.__name__}.{v.__name__}'

    if inspect.isfunction(v):
        return f'{v.__module__}.{v.__qualname__}'

    if hasattr(v,'__class__') and hasattr(v,'__call__'):
        c=v.__class__
        return f'{c.__module__}.{c.__name__}'

    return repr(v)


def _mark_login_not_required(url_patterns):
    for entry in url_patterns:
        if isinstance(entry, URLResolver):
            _mark_login_not_required(entry.url_patterns)
        elif isinstance(entry, URLPattern):
            login_not_required(entry.callback)


def login_not_required_include(urlconf_module: str):
    """ include() for a third-party URLconf (registration, OIDC) whose views must be reachable before login """
    urlconf = import_module(urlconf_module)
    _mark_login_not_required(urlconf.urlpatterns)
    return include(urlconf)
