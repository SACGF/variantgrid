import abc
import logging

from django.core.exceptions import FieldError
from django.core.paginator import Paginator
from django.db import connection

from library.jqgrid.jqgrid_user_row_config import JqGridUserRowConfig
from library.utils import Struct, double_quote
from library.utils.database_utils import iter_dictfetchall


def get_overrides(columns, column_data, model_field=True, queryset_field=True):
    base_colmodel_override = {
        'model_field': True,
        'queryset_field': True,
        'editable': True,
    }
    overrides = {}

    for c, col_data_dict in zip(columns, column_data):
        colmodel = base_colmodel_override.copy()
        colmodel.update(col_data_dict)
        colmodel['name'] = c
        colmodel['index'] = c
        colmodel['search'] = model_field  # Appears in filter dropdown. Currently has to be model due to FK lookup
        colmodel['model_field'] = model_field
        colmodel['queryset_field'] = queryset_field
        overrides[c] = colmodel
    return overrides


class JqGridSQL(JqGridUserRowConfig):
    paginate = True

    @abc.abstractmethod
    def get_count(self, request):
        raise NotImplementedError()

    def is_empty(self, request):
        return self.get_count(request) == 0

    @abc.abstractmethod
    def get_sql_params_and_columns(self, request):
        raise NotImplementedError()

    def get_page_and_paginator(self, request, items):
        paginate_by = self.get_paginate_by(request)
        if not (self.paginate and paginate_by):
            return None, None

        paginator = Paginator(items, paginate_by,
                              allow_empty_first_page=self.allow_empty)

        paginator.count = self.get_count(request)
        page = request.GET.get('page', 1)
        return int(page), paginator

    def get_items(self, request):
        if self.is_empty(request):
            paginator = Struct(num_pages=0, count=0)
            count = Struct(number=0)
            return paginator, count, []
        return self.get_items_from_queryset(request)

    def get_items_from_queryset(self, request):
        sql, params, column_names, is_sorted = self.get_sql_params_and_columns(request)

        if not is_sorted:
            sidx = request.GET.get('sidx')
            if sidx:
                if sidx not in self.fields:
                    fields_comma_sep = ", ".join(sorted(self.fields))
                    msg = f"sidx '{sidx}' is not in fields '{fields_comma_sep}'"
                    raise FieldError(msg)

                # Default to sort by column itself
                override = self.get_override(sidx)
                order_by = override.get("order_by", double_quote(sidx))
                sord = request.GET.get('sord')
                sord_sql = 'DESC' if sord == 'desc' else 'ASC'
                sql += f' ORDER BY {order_by} {sord_sql}'

        page, paginator = self.get_page_and_paginator(request, [])
        if paginator is not None:
            offset = (page - 1) * paginator.per_page
            offset = max(0, offset)  # Don't allow below 0
            sql += f" LIMIT {paginator.per_page} OFFSET {offset}"

        cursor = connection.cursor()
        try:
            # Sometimes we've already escaped the SQL to pull it apart - so will end up having SQL "like '%splice%'"
            # We don't want to interpolate the %s there. In that case params will be empty, so don't pass
            if params:
                cursor.execute(sql, params)
            else:
                cursor.execute(sql)
        except Exception as e:
            logging.error("SQL Error:")
            logging.error("SQL: %s", sql)
            logging.error("params: %s", params)
            raise e

        items = iter_dictfetchall(cursor, column_names)
        items = self.iter_format_items(items)
        page = Struct(number=page)
        return paginator, page, items
