from django.http import HttpRequest
from django.http.response import HttpResponseBase

from classification.views.exports.classification_export_formatter2 import ClassificationExportFormatter2


def serve_export(request: HttpRequest) -> HttpResponseBase:
    type = request.query_params.get('type')
    formatter: ClassificationExportFormatter2
    if type == 'mvl':
        # TODO there has to be a better way to make the formatter automatically?
        from classification.views.exports.classification_export_formatter2_mvl import ClassificationExportFormatter2MVL
        formatter = ClassificationExportFormatter2MVL.from_request(request)
    return formatter.serve()
