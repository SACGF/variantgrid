from typing import Iterable, Optional, Any, Type, Dict, List, Iterator
from dateutil.tz import gettz
from django.conf import settings
from django.http import HttpRequest, StreamingHttpResponse

from library.utils.text_utils import delimited_row
from library.utils.date_utils import local_date_string
import inspect
from datetime import datetime
import json
from library.utils.json_utils import JsonObjType


def export_column(label: Optional[str] = None, sub_data: Optional[Type] = None, categories: Dict[str, Any] = None, format: Dict[str, Any] = None):
    """
    Extend ExportRow and annotate methods with export_column.
    The order of defined methods determines the order that the results will appear in an export file
    :param label: The label that will appear in the CSV header (defaults to method name if not provided)
    :param sub_data: An optional SubType of another ExportRow for nested data
    :param categories: If this export column is only valid in some contexts, provide it here
    :param format: formatters used to idenfify data type, e.g. {"tz": "default"}
    """
    tz_name = settings.TIME_ZONE
    use_tz = None
    if format and (format_tz_name := format.get("tz")):
        if format_tz_name == "default":
            format_tz_name = settings.TIME_ZONE
        use_tz = gettz(format_tz_name)
        tz_name = use_tz.tzname(datetime.now())

    def decorator(method):
        def wrapper(*args, **kwargs):
            result = method(*args, **kwargs)
            if isinstance(result, datetime):
                if use_tz:
                    result = result.astimezone(use_tz)
                return result.strftime("%Y-%m-%d %H:%M")

            return result
        # have to cache the line number of the source method, otherwise we just get the line number of this wrapper
        wrapper.line_number = inspect.getsourcelines(method)[1]
        wrapper.label = label or method.__name__
        if '$site_name' in wrapper.label:
            wrapper.label = wrapper.label.replace('$site_name', settings.SITE_NAME)

        if use_tz:
            wrapper.label = f"{wrapper.label} ({tz_name})"

        wrapper.__name__ = method.__name__
        wrapper.is_export = True
        wrapper.categories = categories
        wrapper.sub_data = sub_data
        return wrapper
    return decorator


class ExportRow:

    @staticmethod
    def get_export_methods(klass, categories: Optional[Dict[str, Any]] = None):
        # TODO, provide export_methods in ExportRow but still work out if we need to search or not
        if not hasattr(klass, 'export_methods'):
            export_methods = [func for _, func in inspect.getmembers(klass, lambda x: getattr(x, 'is_export', False))]
            export_methods.sort(key=lambda x: x.line_number)
            klass.export_methods = export_methods

            if not klass.export_methods:
                raise ValueError(f"ExportRow class {klass} has no @export_columns, make sure you annotate @export_column(), not @export_column")

        export_methods = klass.export_methods
        if categories:
            def passes_filter(export_method) -> bool:
                nonlocal categories
                decorated_values = export_method.categories or dict()
                # for every requirement of categories
                for key, value in categories.items():
                    # get the decorated value
                    if decorated_value := decorated_values.get(key):
                        # handle decorated value being a collection (and matching a single value in that collection)
                        if isinstance(decorated_value, (set, tuple, list)):
                            if value not in decorated_value:
                                return False
                        elif decorated_value != value:
                            return False
                    # if the requirement for the category is None and there's no value at all in the decorator
                    # it passes the test
                    elif value is not None:
                        return False

                return True

            export_methods = [em for em in export_methods if passes_filter(em)]

        return export_methods

    @classmethod
    def _data_generator(cls, data: Iterable[Any]) -> Iterator[Any]:
        for row_data in data:
            if row_data is None:
                continue
            if not isinstance(row_data, cls):
                # it's expected that the class can be initiated with each "row" in data
                row_data = cls(row_data)
            yield row_data

    @classmethod
    def csv_generator(cls, data: Iterable[Any], delimiter=',', include_header=True, categories: Optional[Dict[str, Any]] = None) -> Iterator[str]:
        try:
            if include_header:
                yield delimited_row(cls.csv_header(categories=categories), delimiter=delimiter)
            for row_data in cls._data_generator(data):
                yield delimited_row(row_data.to_csv(categories=categories), delimiter=delimiter)
        except:
            from library.log_utils import report_exc_info
            report_exc_info(extra_data={"activity": "Exporting"})
            yield "** File terminated due to error"
            raise

    @classmethod
    def json_generator(cls, data: Iterable[Any], records_key: str = "records", categories: Optional[Dict[str, Any]] = None) -> Iterator[str]:
        """
        :param data: Iterable data of either cls or that can be passed to cls's constructor
        :param records_key:
        :param categories:
        :return: A generator that will produce several strings, that when concat makes valid Json
        """
        first_row = True
        try:
            if records_key:
                yield f'{{"{records_key}": ['

            for row_data in cls._data_generator(data):
                text = ""
                if not first_row:
                    text += ",\n"
                first_row = False
                text += json.dumps(row_data.to_json(categories=categories))
                yield text

            if records_key:
                yield ']}}'
        except:
            from library.log_utils import report_exc_info
            report_exc_info(extra_data={"activity": "Exporting"})
            yield "\"error\"** File terminated due to error"
            raise

    @classmethod
    def csv_header(cls, categories: Optional[Dict[str, Any]] = None) -> List[str]:
        row = list()
        for method in ExportRow.get_export_methods(cls, categories=categories):
            label = method.label or method.__name__
            if sub_data := method.sub_data:
                sub_header = sub_data.csv_header(categories=categories)
                for sub in sub_header:
                    row.append(label + "." + sub)
            else:
                row.append(label)
        return row

    def to_csv(self, categories: Optional[Dict[str, Any]] = None) -> List[str]:
        row = list()
        for method in ExportRow.get_export_methods(self.__class__, categories=categories):
            result = method(self)
            if sub_data := method.sub_data:
                if result is None:
                    for _ in sub_data.csv_header(categories=categories):
                        row.append("")
                else:
                    row += result.to_csv()
            else:
                row.append(result)
        return row

    def to_json(self, categories: Optional[Dict[str, Any]] = None) -> JsonObjType:
        row = dict()
        for method in ExportRow.get_export_methods(self.__class__, categories=categories):
            result = method(self)
            value: Any
            if result is None:
                value = None
            elif method.sub_data:
                value = result.to_json(categories=categories)
            else:
                value = result

            if value == "":
                value = None
            # row[method.__name__] = value
            row[method.label or method.__name] = value

        return row

    @classmethod
    def streaming(cls, request: HttpRequest, data: Iterable[Any], filename: str, categories: Optional[Dict[str, Any]] = None):
        if request.GET.get('format') == 'json':
            return cls.streaming_json(data, filename, categories=categories)
        else:
            return cls.streaming_csv(data, filename, categories=categories)

    @classmethod
    def streaming_csv(cls, data: Iterable[Any], filename: str, categories: Optional[Dict[str, Any]] = None):
        date_str = local_date_string()

        response = StreamingHttpResponse(cls.csv_generator(data, categories=categories), content_type='text/csv')
        response['Content-Disposition'] = f'attachment; filename="{filename}_{settings.SITE_NAME}_{date_str}.csv"'
        return response

    @classmethod
    def streaming_json(cls, data: Iterable[Any], filename: str, records_key: str = None, categories: Optional[Dict[str, Any]] = None):
        date_str = local_date_string()

        if not records_key:
            records_key = filename.replace(" ", "_")

        response = StreamingHttpResponse(cls.json_generator(data, records_key, categories=categories), content_type='application/json')
        response['Content-Disposition'] = f'attachment; filename="{filename}_{settings.SITE_NAME}_{date_str}.json"'
        return response
