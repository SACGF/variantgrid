from typing import Callable, Dict

from django.http import HttpRequest

from classification.views.exports.classification_export_formatter2 import ClassificationExportFormatter2

ExportFormatterFactory = Callable[[HttpRequest], ClassificationExportFormatter2]


_classification_export_registry: Dict[str, ExportFormatterFactory] = dict()


class UnsupportedExportType(Exception):
    pass


def register_classification_exporter(format_type: str):
    def _wrapper(cls):
        _classification_export_registry[format_type] = cls.from_request
        cls.format_type = format_type
        return cls
    return _wrapper


def classification_exporter_for_request(request: HttpRequest) -> ClassificationExportFormatter2:
    format_type = request.query_params.get('type')
    if factory := _classification_export_registry.get(format_type):
        return factory(request)
    raise UnsupportedExportType(f"No ClassificationExportFormatter2 could be found for '{format_type}'")