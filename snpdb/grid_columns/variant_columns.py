""" Builds the variant grids' columns.

Their column set is per request - it comes from the user's CustomColumnsCollection plus whatever the
analysis node adds - so it is authored as a {field_name: {overrides}} dict merged with
update_dict_of_dict_values, and turned into RichColumns here. The overrides use RichColumn's own
keyword names, plus 'model_field' / 'queryset_field' (which are VariantGridColumn fields).
"""
import logging
from collections.abc import Callable, Iterable, Iterator
from typing import Any, Optional

from django.db.models import fields
from django.utils.timezone import localtime

from library.django_utils import lookup_field_path
from library.log_utils import log_traceback
from snpdb.views.datatable_view import RichColumn

# These grids are wide (dozens of columns) and hold cells with hundreds of links, so the client lays
# them out table-layout: fixed and clips each cell to one line. That needs every column to carry a
# width. @see .variantgrid-datatable in global.scss
DEFAULT_COLUMN_WIDTH = 150

# Django field -> the filter builder's input widget / whether a column summarises quantitatively
_FIELD_DATA_TYPES = [
    (fields.AutoField, "int"),
    (fields.IntegerField, "int"),
    (fields.FloatField, "float"),
    (fields.DateTimeField, "date"),
]

# Keys handled here rather than passed to RichColumn
_NON_COLUMN_OVERRIDES = {"model_field", "queryset_field", "name"}


def column_overrides(columns: list[str], column_data: list[dict],
                     model_field=True, queryset_field=True) -> dict[str, dict]:
    """ {field_name: overrides} for columns that share the same model_field/queryset_field treatment """
    overrides = {}
    for c, col_data_dict in zip(columns, column_data):
        override = {
            "model_field": model_field,
            "queryset_field": queryset_field,
            # Filtering is by Django lookup, so only real model fields can be offered
            "filterable": model_field,
        }
        override.update(col_data_dict)
        overrides[c] = override
    return overrides


def _choices_formatter(choices: dict) -> Callable:
    """ Closure over the loop variable """
    def format_choice(row, field):
        val = row[field]
        return choices.get(val, val)
    return format_choice


def _format_datetime_as_local_time(row, field):
    if value := row[field]:
        value = localtime(value).strftime("%Y-%m-%d %H:%M")
    return value


def _model_field_column_kwargs(model, field_name: str) -> dict[str, Any]:
    """ Label, type and choices taken off the Django field the column selects """
    field = lookup_field_path(model._meta, field_name)
    column_kwargs: dict[str, Any] = {"label": field.verbose_name}

    if isinstance(field, fields.related.ForeignKey):
        # A FK filters against the related object rather than a value - use the full path to the
        # column you actually want
        column_kwargs["filterable"] = False
        return column_kwargs

    for field_type, data_type in _FIELD_DATA_TYPES:
        if isinstance(field, field_type):
            column_kwargs["data_type"] = data_type
            break

    if isinstance(field, fields.BooleanField):
        column_kwargs["filter_choices"] = {"False": "False", "True": "True"}
    elif choices := field.choices:
        column_kwargs["filter_choices"] = dict(choices)
        if isinstance(field, fields.CharField):
            # Show the expanded string rather than the stored code
            column_kwargs["server_side_formatter"] = _choices_formatter(dict(choices))

    if isinstance(field, fields.DateTimeField):
        column_kwargs["server_side_formatter"] = _format_datetime_as_local_time
    return column_kwargs


def build_variant_columns(model, field_names: list[str], overrides: dict[str, dict]) -> list[RichColumn]:
    """ One RichColumn per field, in field order - hidden columns included, as the client's column
        index has to map back to the column it came from """
    columns = []
    for field_name in field_names:
        override = overrides.get(field_name, {})
        try:
            column_kwargs: dict[str, Any] = {}
            if override.get("model_field", True):
                column_kwargs = _model_field_column_kwargs(model, field_name)
            column_kwargs.update({k: v for k, v in override.items() if k not in _NON_COLUMN_OVERRIDES})
            column_kwargs.setdefault("orderable", True)
            column_kwargs.setdefault("width", DEFAULT_COLUMN_WIDTH)
            column_kwargs.setdefault("filterable", True)
            columns.append(RichColumn(key=field_name, name=field_name,
                                      queryset_field=override.get("queryset_field", True),
                                      include_in_csv=True, search=False,
                                      # Rows with no value for the column sort as the lowest
                                      sort_nulls_first=True, **column_kwargs))
        except Exception:
            logging.error("Field_name: '%s'", field_name)
            log_traceback()
            raise
    return columns


def column_formatters(columns: Iterable[RichColumn]) -> dict[str, Callable]:
    return {rc.name: rc.server_side_formatter for rc in columns if rc.server_side_formatter}


def format_rows(columns: Iterable[RichColumn], items: Iterable[dict]) -> Iterator[dict]:
    """ Applies the columns' server side formatters to each row, in place. Shared by the grid and the
        CSV/VCF export, so the two can't disagree about a value """
    formatters = column_formatters(columns)
    if not formatters:
        yield from items
        return

    rows = iter(items)
    try:
        first_row = next(rows)
    except StopIteration:
        return

    # Rows are .values() dicts so they all have the same keys - try each formatter against the first
    # row to see which ones apply, rather than raising and swallowing a KeyError per absent field on
    # every row. Bound per call as callers pass different row shapes.
    applicable = []
    for field, formatter in formatters.items():
        try:
            first_row[field] = formatter(first_row, field)
        except KeyError:
            continue  # field may not be in columns returned...
        except ValueError:
            pass  # formatter applies, it just can't handle this value
        applicable.append((field, formatter))
    yield first_row

    for row in rows:
        for field, formatter in applicable:
            try:
                row[field] = formatter(row, field)
            except ValueError:
                pass  # formatter can't handle this value
        yield row


def column_by_name(columns: Iterable[RichColumn], name: str) -> Optional[RichColumn]:
    for rc in columns:
        if rc.name == name:
            return rc
    return None
