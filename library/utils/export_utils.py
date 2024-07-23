import inspect
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
from itertools import chain
from typing import Iterable, Optional, Any, Type, Iterator, Callable, Protocol

from dateutil.tz import gettz
from django.conf import settings
from django.http import HttpRequest, StreamingHttpResponse

from library.utils.date_utils import local_date_string
from library.utils.json_utils import JsonObjType
from library.utils.text_utils import delimited_row


class ExportDataType(str, Enum):
    datetime_notz = "datetime_notz"  # if timezone doesn't really matter
    datetime = "datetime"
    date = "date"
    any = "any"


class ExportCellMethod(Protocol):
    line_number: int
    label: str
    sub_data: Optional[Type['ExportRow']]
    categories: Optional[dict[Any, Any]]
    data_type: ExportDataType
    def __call__(self, *args) -> Any: ...


def export_column(
        label: Optional[str] = None,
        sub_data: Optional[Type] = None,
        categories: dict[str, Any] = None,
        data_type: ExportDataType = ExportDataType.any):
    """
    Extend ExportRow and annotate methods with export_column.
    The order of defined methods determines the order that the results will appear in an export file
    :param label: The label that will appear in the CSV header (defaults to method name if not provided)
    :param sub_data: An optional SubType of another ExportRow for nested data
    :param categories: If this export column is only valid in some contexts, provide it here
    :param data_type: Type of data to be printed
    """

    def decorator(method):
        def wrapper(*args, **kwargs):
            return method(*args, **kwargs)

        # have to cache the line number of the source method, otherwise we just get the line number of this wrapper
        wrapper.line_number = inspect.getsourcelines(method)[1]
        wrapper.label = label or method.__name__
        if '$site_name' in wrapper.label:
            wrapper.label = wrapper.label.replace('$site_name', settings.SITE_NAME)

        wrapper.data_type = data_type
        wrapper.__name__ = method.__name__
        wrapper.is_export = True
        wrapper.categories = categories
        wrapper.sub_data = sub_data
        return wrapper
    return decorator


class ExportFormat(str, Enum):
    csv = "csv"
    json = "json"


class ExportSettings:

    RE_TZ_FORMAT = re.compile(r".*/(.*/.*)'[)]")
    UTC = gettz("UTC")

    def __init__(self, tz: timezone):
        self.tz = tz

    def format_heading(self, export_format: ExportFormat, data_type: ExportDataType, label: str) -> str:
        if export_format == ExportFormat.csv:
            if data_type == ExportDataType.datetime:
                tz_name = str(self.tz)
                if match := ExportSettings.RE_TZ_FORMAT.match(tz_name):
                    tz_name = match.group(1)

                return f"{label} ({tz_name})"
        return label

    def format_value(self, export_format: ExportFormat, data_type: ExportDataType, value: Any) -> Any:
        if data_type == ExportDataType.date and isinstance(value, datetime):
            return value.strftime('%Y-%m-%d')
        if export_format == ExportFormat.csv:
            if data_type == ExportDataType.datetime and isinstance(value, datetime):
                value = value.astimezone(self.tz)
                # want to put the timezone in the string %z or %Z BUT... Excel doesn't support that
                # because the one constant through history is Excel ruins everything
                return value.strftime("%Y-%m-%d %H:%M")
        elif export_format == ExportFormat.json:
            # convert to UTF
            if isinstance(value, datetime):
                value = value.astimezone(self.tz)
                return value.isoformat()

        return value

    @staticmethod
    def get_for_request():
        from snpdb.user_settings_manager import UserSettingsManager
        return ExportSettings(tz=UserSettingsManager.get_user_timezone())


@dataclass
class ExportTweak:
    """
    Provides some high level changes (primarily to CSV export)
    """
    sub_label_joiner: Callable[[str, str], str] = lambda a, b: a + "." + b
    """
    A method that combine labels and sub labels
    """

    categories: Optional[dict[Any, Any]] = None
    """
    key, value pairs that will determine which columns are included
    """

    force_method_heading: bool = False
    """
    If true, use method names rather than export column labels (more used internally for making mappings)
    """


ExportTweak.DEFAULT = ExportTweak()


class _ColumnZipperRow:
    """
    Class that helps us combine columns with sub data together (if desired)
    Used for headers and regular rows
    """

    def __init__(self, sub_data_zip_types: set[type]):
        self.row = []
        self.sub_data_zip_types = sub_data_zip_types
        self.sub_data_zips: list[list] = []
        self.last_sub_data_zip_type: Optional[type] = None

    def _apply_sub_data_zips(self):
        if self.sub_data_zips:
            flatten = list(chain(*zip(*self.sub_data_zips)))
            self.row.extend(flatten)
        self.last_sub_data_zip_type = None
        self.sub_data_zips.clear()

    def add_cell(self, value):
        self._apply_sub_data_zips()
        self.row.append(value)

    def add_sub_data(self, sub_data_type: type, values: list):
        if sub_data_type in self.sub_data_zip_types:
            if self.last_sub_data_zip_type is not None and sub_data_type != self.last_sub_data_zip_type:
                # more zipped sub data, but different type
                self._apply_sub_data_zips()

            self.sub_data_zips.append(values)
        else:
            self._apply_sub_data_zips()
            self.row.extend(values)

    def get_row(self) -> list:
        self._apply_sub_data_zips()
        return self.row


def get_decorated_methods(cls, categories: Optional[dict[Any, Any]], attribute: str) -> list[Callable]:
    if not hasattr(cls, 'export_methods'):
        export_methods = [func for _, func in inspect.getmembers(cls, lambda x: getattr(x, attribute, False))]
        export_methods.sort(key=lambda x: x.line_number)
        cls.export_methods = export_methods

        if not cls.export_methods:
            raise ValueError(
                f"Decorated class {cls} has no methods with {attribute}, did you decorate the methods correctly?")

    export_methods = cls.export_methods
    if categories:
        def passes_filter(export_method) -> bool:
            nonlocal categories
            # FIXME, some non-NONE but falsey values could get confused here, e.g. 0 and False
            decorated_values = export_method.categories or {}
            # for every requirement of categories
            for key, value in categories.items():
                # get the decorated value
                value_set: set
                if isinstance(value, (set, tuple, list)):
                    value_set = set(value)
                else:
                    value_set = {value}

                if decorated_value := decorated_values.get(key):
                    decorated_set: set
                    if isinstance(decorated_value, (set, tuple, list)):
                        decorated_set = set(decorated_value)
                    else:
                        decorated_set = {decorated_value}
                    return bool(value_set.intersection(decorated_set))

                    # handle decorated value being a collection (and matching a single value in that collection)
                # if the requirement for the category is None and there's no value at all in the decorator
                # it passes the test
                elif value is not None:
                    return False

            return True

        export_methods = [em for em in export_methods if passes_filter(em)]

    return export_methods


class ExportRow:

    @classmethod
    def get_export_methods(cls, export_tweak: ExportTweak = ExportTweak.DEFAULT) -> list[ExportCellMethod]:
        return get_decorated_methods(cls, categories=export_tweak.categories, attribute="is_export")

    @classmethod
    def _data_generator(cls: Type, data: Iterable[Any]) -> Iterator[Any]:
        for row_data in data:
            if row_data is None:
                continue
            if not isinstance(row_data, cls):
                # it's expected that the class can be initiated with each "row" in data
                row_data = cls(row_data)
            yield row_data

    @classmethod
    def csv_generator(cls, data: Iterable[Any], delimiter=',', include_header=True, export_tweak: ExportTweak = ExportTweak.DEFAULT, export_settings: Optional[ExportSettings] = None) -> Iterator[str]:
        if not export_settings:
            export_settings = ExportSettings.get_for_request()
        try:
            if include_header:
                yield delimited_row(cls.csv_header(export_tweak=export_tweak, export_settings=export_settings), delimiter=delimiter)
            for row_data in cls._data_generator(data):
                yield delimited_row(row_data.to_csv(export_tweak=export_tweak, export_settings=export_settings), delimiter=delimiter)
        except:
            from library.log_utils import report_exc_info
            report_exc_info(extra_data={"activity": "Exporting"})
            yield "** File terminated due to error"
            raise

    # Warning: We don't actually do anything with the formatting of JSON objects

    @classmethod
    def json_generator(cls, data: Iterable[Any], records_key: str = "records", export_tweak: ExportTweak = ExportTweak.DEFAULT) -> Iterator[str]:
        """
        :param data: Iterable data of either cls or that can be passed to cls's constructor
        :param records_key:
        :param export_tweak:
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
                text += json.dumps(row_data.to_json(export_tweak=export_tweak))
                yield text

            if records_key:
                yield ']}}'
        except:
            from library.log_utils import report_exc_info
            report_exc_info(extra_data={"activity": "Exporting"})
            yield "\"error\"** File terminated due to error"
            raise

    @classmethod
    def zip_sub_data(cls) -> set[type]:
        return set()

    @classmethod
    def csv_header(cls, export_tweak: ExportTweak = ExportTweak.DEFAULT, export_settings: Optional[ExportSettings] = None) -> list[str]:
        # a bit of duplication between generating the header
        # and generating the data
        # consider making a model where we calculate cells (in the correct order) and can then query those cells
        # for labels or values (given an input)
        if not export_settings:
            export_settings = ExportSettings.get_for_request()

        row = _ColumnZipperRow(sub_data_zip_types=cls.zip_sub_data())

        for method in cls.get_export_methods(export_tweak=export_tweak):
            label: str
            if export_tweak.force_method_heading:
                label = method.__name__
            else:
                label = method.label or method.__name__
                label = export_settings.format_heading(export_format=ExportFormat.csv, data_type=method.data_type, label=label)

            if sub_data := method.sub_data:
                sub_header = sub_data.csv_header(export_tweak=export_tweak, export_settings=export_settings)
                sub_headings = [export_tweak.sub_label_joiner(label, sub) for sub in sub_header]
                row.add_sub_data(sub_data, sub_headings)
            else:
                row.add_cell(label)
        return row.get_row()

    def to_csv(self, export_tweak: ExportTweak = ExportTweak.DEFAULT, export_settings: Optional[ExportSettings] = None) -> list[str]:
        if not export_settings:
            export_settings = ExportSettings.get_for_request()

        row = _ColumnZipperRow(sub_data_zip_types=self.zip_sub_data())

        for method in self.__class__.get_export_methods(export_tweak=export_tweak):
            result = method(self)
            if sub_data := method.sub_data:
                sub_data_values: list
                if result is None:
                    sub_data_values = [""] * len(sub_data.csv_header(export_tweak=export_tweak, export_settings=export_settings))
                else:
                    sub_data_values = result.to_csv(export_tweak=export_tweak, export_settings=export_settings)

                row.add_sub_data(sub_data, sub_data_values)
            else:
                result = export_settings.format_value(export_format=ExportFormat.csv, data_type=method.data_type, value=result)
                row.add_cell(result)
        return row.get_row()

    def to_json(self, export_tweak: ExportTweak = ExportTweak.DEFAULT, export_settings: Optional[ExportSettings] = None) -> JsonObjType:
        if not export_settings:
            export_settings = ExportSettings.get_for_request()
        row = {}
        for method in self.__class__.get_export_methods(export_tweak=export_tweak):
            result = method(self)
            value: Any
            if result is None:
                value = None
            elif method.sub_data:
                value = result.to_json(export_tweak=export_tweak, export_settings=export_settings)
            else:
                value = result

            if value == "":
                value = None

            value = export_settings.format_value(export_format=ExportFormat.json, data_type=method.data_type, value=value)
            row[method.label] = value

        return row

    @classmethod
    def streaming(cls, request: HttpRequest, data: Iterable[Any], filename: str, export_tweak: ExportTweak = ExportTweak.DEFAULT):
        if request.GET.get('format') == 'json':
            return cls.streaming_json(data, filename, export_tweak=export_tweak)
        else:
            return cls.streaming_csv(data, filename, export_tweak=export_tweak)

    @classmethod
    def streaming_csv(cls, data: Iterable[Any], filename: str, export_tweak: ExportTweak = ExportTweak.DEFAULT):
        date_str = local_date_string()

        response = StreamingHttpResponse(cls.csv_generator(data, export_tweak=export_tweak), content_type='text/csv')
        response['Content-Disposition'] = f'attachment; filename="{filename}_{settings.SITE_NAME}_{date_str}.csv"'
        return response

    @classmethod
    def streaming_json(cls, data: Iterable[Any], filename: str, records_key: str = None, export_tweak: ExportTweak = ExportTweak.DEFAULT):
        date_str = local_date_string()

        if not records_key:
            records_key = filename.replace(" ", "_")

        response = StreamingHttpResponse(cls.json_generator(data, records_key, export_tweak=export_tweak), content_type='application/json')
        response['Content-Disposition'] = f'attachment; filename="{filename}_{settings.SITE_NAME}_{date_str}.json"'
        return response
