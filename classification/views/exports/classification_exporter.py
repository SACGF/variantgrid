from django.http import HttpRequest
from django.http.response import HttpResponseBase

from classification.views.exports.classification_export_formatter2 import ClassificationExportFormatter2, \
    ClassificationFilter


def serve_export(request: HttpRequest) -> HttpResponseBase:
    type = request.query_params.get('type')
    filter = ClassificationFilter.from_request(request)
    formatter: ClassificationExportFormatter2
    if type == 'mvl':
        filter.row_limit = 9999

        from classification.views.exports.classification_export_formatter2_mvl import MVLFormatDetails, \
            ClassificationExportFormatter2MVL

        format = MVLFormatDetails.from_request(request)
        formatter = ClassificationExportFormatter2MVL(filter, format)
    return formatter.serve()
