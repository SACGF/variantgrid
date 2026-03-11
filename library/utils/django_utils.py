import inspect
from typing import Optional, TypeVar

from cache_memoize import cache_memoize
from django.conf import settings
from django.contrib.contenttypes.models import ContentType
from django.db.models import Model
from django.http import HttpRequest, HttpResponse
from django.shortcuts import render
from django.template.loader import render_to_string

from library.git import Git


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


ModelT = TypeVar("ModelT", bound=Model)


def refresh_for_update(obj: ModelT) -> ModelT:
    return type(obj).objects.filter(pk=obj.pk).select_for_update().get()


def get_model_content_type_dict(model):
    content_type = ContentType.objects.get_for_model(model)
    return {
        'app_label': content_type.app_label,
        'model': content_type.model
    }


@cache_memoize(30)
def get_cached_project_git_hash() -> str:
    return Git(settings.BASE_DIR).hash


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