from typing import Dict, Optional

from django.http import HttpRequest, HttpResponse
from django.shortcuts import render
from django.template.loader import render_to_string


def render_ajax_view(
        request: HttpRequest,
        template_name: str,
        context: Dict,
        menubar: Optional[str] = None) -> HttpResponse:
    if not context:
        context = {}
    if request.META.get('HTTP_X_REQUESTED_WITH') == 'XMLHttpRequest':
        context['render_mode'] = 'ajax'
        return render(request, template_name, context)
    else:
        context['render_mode'] = 'embedded'
        text = render_to_string(template_name, context, request=request)
        return render(request, "snpdb/embedded_ajax.html", {"embedded_content": text, "menubar": menubar})
