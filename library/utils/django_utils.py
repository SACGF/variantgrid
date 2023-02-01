from typing import Dict

from django.http import HttpRequest, HttpResponse
from django.shortcuts import render
from django.template.loader import render_to_string


def render_ajax_view(request: HttpRequest, template_name: str, context: Dict) -> HttpResponse:
    if request.META.get('HTTP_X_REQUESTED_WITH') == 'XMLHttpRequest':
        return render(request, template_name, context)
    else:
        text = render_to_string(template_name, context, request=request)
        return render(request, "snpdb/embedded_ajax.html", {"embedded_content": text})
