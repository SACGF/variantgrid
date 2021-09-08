# Copyright (c) 2009, Gerry Eisenhaur
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#    1. Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#
#    2. Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#
#    3. Neither the name of the project nor the names of its contributors may
#       be used to endorse or promote products derived from this software
#       without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import json
import logging
import operator
from copy import deepcopy
from functools import reduce

from django.core.exceptions import FieldError, ImproperlyConfigured
from django.core.paginator import Paginator, InvalidPage
from django.core.serializers.json import DjangoJSONEncoder
from django.db.models import fields, JSONField, F
from django.db.models.fields.json import KeyTransform
from django.db.models.query_utils import Q
from django.utils.encoding import smart_str
from django.utils.timezone import localtime

from library.log_utils import log_traceback


def json_encode(data):
    encoder = DjangoJSONEncoder(allow_nan=False)  # strict JSON so JS doesn't crash
    return encoder.encode(data)


def format_operation(op):
    ops = {
        "eq": "equal",
        "ne": "not equal",
        "bw": "begins with",
        "bn": "does not begin with",
        "ew": "ends with",
        "en": "does not end with",
        "cn": "contains",
        "nc": "does not contain",
        "nu": "is null",
        "nn": "is not null",
        "in": "is in",
        "ni": "is not in",
        'gt': "greater than",
        'ge': "greater than or equal to",
        'lt': "less than",
        'le': "less than or equal to",
    }
    return ops[op]


class JqGrid:
    queryset = None
    model = None
    fields = []
    allow_empty = True

    pager_id = '#pager'
    url = None
    caption = None
    colmodel_overrides = {}  # Class member, do not modify in instances!

    def __init__(self):
        self.extra_config = {}
        self._overrides = deepcopy(self.colmodel_overrides)
        if not self.fields:
            msg = "You need to populate fields before calling __init__"
            raise ValueError(msg)

        # Colmodel indexes for foreign keys (allows sorting by FK)
        for f in self.fields:
            if f.find('__') > 0:
                override = self._overrides.get(f, {})
                override['index'] = f
                self._overrides[f] = override

    def update_overrides(self, override: dict):
        for k, v in override.items():
            old_override = self._overrides.get(k, {})
            old_override.update(v)
            self._overrides[k] = old_override

    def get_override(self, field_name):
        return self._overrides.get(field_name, {})

    def add_annotations_for_json_fields(self, queryset):
        # You can't join into a JSON field like values("json_field__key")
        # normally, but can use a trick see https://stackoverflow.com/a/45369944/295724

        # You annotate into the query, field by field, eg "json_field__json_key_0__json_key_1", ops will be:
        # qs = qs.annotate(json_key_0__json_key_1=KeyTransform("json_key_0", "json_field"))
        # qs = qs.annotate(json_field__json_key_0__json_key_1=KeyTransform("json_key_1", "json_key_0__json_key_1")

        opts = self.get_model()._meta
        for f in self.fields:
            data = self.lookup_foreign_key_field(opts, f, return_json_fields_as_tuple=True)
            if isinstance(data, tuple):
                (_, json_data_field) = data
                ri = f.rindex("__" + json_data_field)
                json_field = f[:ri]  # strip json data fields to get json field
                json_data_field_list = json_data_field.split("__")
                full_field_list = [json_field] + json_data_field_list
                key_field = json_field
                i = 0
                num_json_fields = len(json_data_field_list)
                while i < num_json_fields:
                    value_field = json_data_field_list[i]
                    back_index = num_json_fields - i
                    anno_field = '__'.join(full_field_list[back_index - 1:])
                    kwargs = {anno_field: KeyTransform(value_field, key_field)}
                    #print("i=%d: kwargs=%s" % (i, kwargs))
                    queryset = queryset.annotate(**kwargs)
                    key_field = anno_field
                    i += 1

        return queryset

    def get_queryset(self, request):
        if hasattr(self, 'queryset') and self.queryset is not None:
            queryset = self.queryset._clone()
        elif hasattr(self, 'model') and self.model is not None:
            queryset = self.model.objects.all()
            queryset = self.add_annotations_for_json_fields(queryset)
            queryset = queryset.values(*self.get_field_names())
        else:
            raise ImproperlyConfigured("No queryset or model defined.")
        self.queryset = queryset
        return self.queryset

    def get_model(self):
        if hasattr(self, 'model') and self.model is not None:
            model = self.model
        elif hasattr(self, 'queryset') and self.queryset is not None:
            model = self.queryset.model
            self.model = model
        else:
            raise ImproperlyConfigured("No queryset or model defined.")
        return model

    def get_items(self, request):
        items = self.get_queryset(request)
        items = self.filter_items(request, items)
        items = self.sort_items(request, items)
        paginator, page, items = self.paginate_items(request, items)
        items = self.iter_format_items(items)
        return paginator, page, items

    def get_filters(self, request):
        _search = request.GET.get('_search')
        filters = None

        if _search == 'true':
            _filters = request.GET.get('filters')
            try:
                filters = _filters and json.loads(_filters)
            except ValueError:
                return None

            if filters is None:
                field = request.GET.get('searchField')
                op = request.GET.get('searchOper')
                data = request.GET.get('searchString')

                if all([field, op, data]):
                    filters = {
                        'groupOp': 'AND',
                        'rules': [{'op': op, 'field': field, 'data': data}]
                    }
        return filters

    def filter_items(self, request, items):
        _filters = self.get_filters(request)
        if _filters:
            q = self.get_q(_filters)
            if q:
                return items.filter(q)

        return items

    def get_q(self, json_filters):
        # TODO: Add more support for RelatedFields (searching and displaying)
        # FIXME: Validate data types are correct for field being searched.
        filter_map = {
            # jqgrid op: (django_lookup, use_exclude)
            'ne': ('%(field)s__exact', True),
            'bn': ('%(field)s__startswith', True),
            'en': ('%(field)s__endswith', True),
            'nc': ('%(field)s__contains', True),
            'ni': ('%(field)s__in', True),
            'in': ('%(field)s__in', False),
            'eq': ('%(field)s__exact', False),
            'bw': ('%(field)s__startswith', False),
            'gt': ('%(field)s__gt', False),
            'ge': ('%(field)s__gte', False),
            'lt': ('%(field)s__lt', False),
            'le': ('%(field)s__lte', False),
            'ew': ('%(field)s__endswith', False),
            'cn': ('%(field)s__contains', False),
            'nu': ('%(field)s__isnull', True),
            'nn': ('%(field)s__isnull', False),
        }
        if self.get_config(False)['ignoreCase']:
            filter_map.update({'ne': ('%(field)s__iexact', True),
                               'eq': ('%(field)s__iexact', False),
                               'bn': ('%(field)s__istartswith', True),
                               'bw': ('%(field)s__istartswith', False),
                               'en': ('%(field)s__iendswith', True),
                               'ew': ('%(field)s__iendswith', False),
                               'nc': ('%(field)s__icontains', True),
                               'cn': ('%(field)s__icontains', False)
                               }
                              )

        q_filters = []
        for rule in json_filters['rules']:
            op, field, data = rule['op'], rule['field'], rule['data']

            filter_fmt, exclude = filter_map[op]
            filter_str = smart_str(filter_fmt % {'field': field})
            if filter_fmt.endswith('__in'):
                filter_kwargs = {filter_str: data.split(',')}
            elif filter_fmt.endswith('__isnull'):
                # FilterNode was slow - pass exclude as arg and don't invert Q
                # Prev generated code like:
                # (NOT (AND: ('variantannotation__dbsnp_rs_id__isnull', True))
                # which generated a full table scan on variant (>100M+ rows...)
                filter_kwargs = {filter_str: exclude}
                exclude = False
            else:
                filter_kwargs = {filter_str: smart_str(data)}

            q = Q(**filter_kwargs)
            if exclude:
                q = ~q
            q_filters.append(q)

        filters = []
        if q_filters:  # Check for search w/no rules
            if json_filters['groupOp'].upper() == 'OR':
                filters = reduce(operator.or_, q_filters)
            else:
                filters = reduce(operator.and_, q_filters)
        return filters

    def sort_items(self, request, items):
        sidx = request.GET.get('sidx')
        if sidx is not None:
            sord = request.GET.get('sord')
            order_by = F(sidx)
            # dlawrence - sort nulls first/last via
            # https://docs.djangoproject.com/en/3.1/ref/models/expressions/#using-f-to-sort-null-values
            if sord == "desc":
                order_by = order_by.desc(nulls_last=True)
            else:
                order_by = order_by.asc(nulls_first=True)

            # always fall back to a second sort column
            # to ensure reliable results, by default have that
            # column be descending id (so newest are first)

            second_sidx = 'pk'
            second_order_by = '-pk'

            if 'sortname' in self.extra_config:
                second_sidx = self.extra_config.get('sortname')
                second_sort_order = ('-' if self.extra_config.get('sortorder') == 'desc' else '')
                second_order_by = f"{second_sort_order}{second_sidx}"

            try:
                if sidx == second_sidx:
                    items = items.order_by(order_by)
                else:
                    items = items.order_by(order_by, second_order_by)
            except FieldError as fe:
                print(fe)
        return items

    def get_paginate_by(self, request):
        rows = request.GET.get('rows', self.get_config(False)['rowNum'])
        try:
            paginate_by = int(rows)
        except ValueError:
            paginate_by = 10
        return paginate_by

    def paginate_items(self, request, items):
        paginate_by = self.get_paginate_by(request)
        if not paginate_by:
            return None, None, items

        paginator = Paginator(items, paginate_by,
                              allow_empty_first_page=self.allow_empty)
        page = request.GET.get('page', 1)

        try:
            page_number = int(page)
            page = paginator.page(page_number)
        except (ValueError, InvalidPage):
            page = paginator.page(1)
        return paginator, page, page.object_list

    def iter_format_items(self, items):
        """ Replace CharField choices with expanded strings
            Dates with formatted local time
        """

        field_formatters = self.get_field_formatters()
        if field_formatters:

            def iter_formatted_items(rows):
                for row in rows:
                    for f, formatter in field_formatters.items():
                        try:
                            row[f] = formatter(row, f)
                        except (KeyError, ValueError):
                            pass  # field may not be in columns returned...
                    yield row

            items = iter_formatted_items(items)
        return items

    def get_data(self, request):
        paginator, page, items = self.get_items(request)
        return {
            'page': int(page.number),
            'total': int(paginator.num_pages),
            'rows': list(items),
            'records': int(paginator.count),
        }

    def get_json(self, request):
        data = self.get_data(request)
        return json_encode(data)

    def get_default_config(self):
        config = {
            'datatype': 'json',
            'autowidth': True,
            'forcefit': True,
            'ignoreCase': True,
            'shrinkToFit': True,
            'jsonReader': {'repeatitems': False},
            'rowNum': 10,
            'rowList': [10, 25, 50, 100],
            'sortname': 'pk',
            'viewrecords': True,
            'sortorder': "asc",
            'pager': self.pager_id,
            'altRows': True,
            'gridview': True,
            'height': 'auto',
            #'multikey': 'ctrlKey',
            #'multiboxonly': True,
            #'multiselect': True,
            #'toolbar': [False, 'bottom'],
            #'userData': None,
            #'rownumbers': False,
        }
        return config

    def get_url(self):
        return str(self.url)

    def get_caption(self):
        if self.caption is None:
            opts = self.get_model()._meta
            self.caption = opts.verbose_name_plural.capitalize()
        return self.caption

    def get_config(self, as_json=True):
        config = self.get_default_config()
        config.update(self.extra_config)
        config.update({
            'url': self.get_url(),
            'caption': self.get_caption(),
            'colModel': self.get_colmodels(remove_server_side_only=True),
        })
        if as_json:
            config = json_encode(config)
        return config

    @staticmethod
    def lookup_foreign_key_field(options, field_name, return_json_fields_as_tuple=False):
        """ Make a field lookup converting __ into real models fields
            returns field normally.
            If return_json_fields_as_tuple=True, then return strings of (json_field, json_data_field)
        """

        if '__' in field_name:
            fk_name, field_name = field_name.split('__', 1)
            # dlawrence: Allows you to go into one to one model relationships
            field = options.get_field(fk_name)
            if field:
                if isinstance(field, JSONField):
                    # Don't descend into JSON fields
                    if return_json_fields_as_tuple:
                        return fk_name, field_name
                    return field
                foreign_model_options = field.related_model._meta
            return JqGrid.lookup_foreign_key_field(foreign_model_options, field_name,
                                                   return_json_fields_as_tuple=return_json_fields_as_tuple)
        return options.get_field(field_name)

    def get_colmodels(self, remove_server_side_only=False):
        colmodels = []
        opts = self.get_model()._meta
        for field_name in self.get_field_names():
            try:
                override = self.get_override(field_name)
                if override.get("model_field", True):
                    field = self.lookup_foreign_key_field(opts, field_name)
                    colmodel = JqGrid.field_to_colmodel(field, field_name)
                else:
                    colmodel = {}

                if override:
                    if remove_server_side_only:
                        override = override.copy()
                        override.pop('server_side_formatter', None)
                    colmodel.update(override)

                colmodels.append(colmodel)
            except:
                logging.error("Field_name: '%s'", field_name)
                log_traceback()
                raise
        return colmodels

    def get_field_choice(self, field, field_name):
        field_choice = None
        if isinstance(field, fields.CharField):
            if field.choices:
                field_choice = dict(field.choices)
        return field_choice

    def get_field_choices(self):
        field_choices = {}
        opts = self.get_model()._meta
        for field_name in self.get_field_names():
            override = self.get_override(field_name)
            if override.get("model_field", True):
                field = self.lookup_foreign_key_field(opts, field_name)
                field_choice = self.get_field_choice(field, field_name)
            else:
                field_choice = {}

            if field_choice:
                field_choices[field_name] = field_choice
        return field_choices

    def get_datetime_fields(self):
        datetime_fields = []
        opts = self.get_model()._meta
        for field_name in self.get_field_names():
            override = self.get_override(field_name)
            if override.get("model_field", True):
                field = self.lookup_foreign_key_field(opts, field_name)
                if isinstance(field, fields.DateTimeField):
                    datetime_fields.append(field_name)

        return datetime_fields

    @staticmethod
    def _make_choices_formatter(choices):
        """ Need this to perform a closure over a loop variable  """
        def choices_formatter(row, field):
            val = row[field]
            return choices.get(val, val)
        return choices_formatter

    def get_field_formatters(self):
        """ returns dict of column names with a formatter function """
        field_formatters = {}
        for field_name, choices in self.get_field_choices().items():
            override = self.get_override(field_name)
            choice_display = override.get("choice_display", True)  # On by default
            if choice_display:
                field_formatters[field_name] = self._make_choices_formatter(choices)

        def date_format_as_local_time(row, field):
            value = row[field]
            if value:
                lt = localtime(value)
                value = lt.strftime("%Y-%m-%d %H:%M")
            return value

        for c in self.get_datetime_fields():
            field_formatters[c] = date_format_as_local_time

        for field_name in self.get_field_names():
            override = self.get_override(field_name)
            server_side_formatter = override.get("server_side_formatter")
            if server_side_formatter:
                field_formatters[field_name] = server_side_formatter

        return field_formatters

    def get_field_names(self):
        """ Returns self.fields if set, otherwise model fields """
        fields = list(self.fields)  # So caller won't modify internals
        if not fields:
            fields = [f.name for f in self.get_model()._meta.local_fields]
        return fields

    @staticmethod
    def add_search_options(field, colmodel):
        # See: http://www.trirand.com/jqgridwiki/doku.php?id=wiki:search_config
        is_boolean_field = any([isinstance(field, t) for t in [fields.BooleanField, fields.NullBooleanField]])
        if is_boolean_field:
            colmodel['stype'] = 'select'
            colmodel['searchoptions'] = {'value': {'False': 'False', 'True': 'True'}}

        if field.choices:
            colmodel['stype'] = 'select'
            colmodel['searchoptions'] = {'value': dict(field.choices)}

    @staticmethod
    def field_to_colmodel(field, field_name):
        colmodel_rules = {"required": not (field.null and field.blank)}

        colmodel = {
            'name': field_name,
            'index': field_name,  # Index is server side column name
            'label': field.verbose_name,
            'editable': True,
            'editrules': colmodel_rules,
            'searchrules': colmodel_rules,
        }

        # VariantGrid #182 - Searching rules to stop type errors
        types = [(fields.AutoField, {'sorttype': 'int', 'rule': 'integer'}),
                 (fields.IntegerField, {'sorttype': 'int', 'rule': 'integer'}),
                 (fields.FloatField, {'sorttype': 'float', 'rule': 'number'}),
                 (fields.DateTimeField, {'sorttype': 'date', 'rule': 'date'})]

        if isinstance(field, fields.related.ForeignKey):
            logging.warning("JQGrid colmodel field '%s' is ForeignKey - disabling search. "
                            "Use full path to final column to avoid filter errors", field)
            colmodel["search"] = False
        else:
            for t, type_info in types:
                if isinstance(field, t):
                    sorttype = type_info.get('sorttype')
                    if sorttype:
                        colmodel['sorttype'] = sorttype

                    rule = type_info.get('rule')
                    if rule:
                        colmodel_rules[rule] = True
                    break

            JqGrid.add_search_options(field, colmodel)

        return colmodel
