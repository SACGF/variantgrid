from django.http import HttpRequest
from django.http.response import HttpResponseBase
from classification.views.exports.classification_export_decorator import classification_exporter_for_request


def serve_export(request: HttpRequest) -> HttpResponseBase:
    return classification_exporter_for_request(request).serve()