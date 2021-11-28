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

import abc
import json
import logging
import math
from copy import deepcopy

# TODO: Merge code from jqGrid and this together?
from library.pandas_utils import df_nan_to_none


class DataFrameJqGrid:
    fields = []
    allow_empty = True

    pager_id = '#pager'
    url = None
    caption = None
    colmodel_overrides = {}

    @abc.abstractmethod
    def get_dataframe(self):
        pass

    def __init__(self):
        self.extra_config = {}
        self._overrides = deepcopy(self.colmodel_overrides)

    def get_override(self, field_name):
        return self._overrides.get(field_name, {})

    def get_items(self, request):
        df = self.get_dataframe()
        df = df_nan_to_none(df)  # JSON can't handle NaN
        df = self.sort_dataframe(request, df)
        (df_slice, page_number, num_pages, num_records) = self.paginate_df(request, df)
        items = self.get_dataframe_rows(df_slice)
        return page_number, num_pages, items, num_records

    def sort_dataframe(self, request, df):
        sidx = request.GET.get('sidx')
        if sidx is not None:
            ascending = not request.GET.get('sord') == 'desc'
            if sidx == 'ID':
                df = df.sort_index(ascending=ascending)
            elif sidx in df.columns:
                df = df.sort_values(sidx, ascending=ascending)
            else:
                logging.error("Request: %s asked for grid sort index (sidx) of '%s' which is not in DF!", request.GET, sidx)
        return df

    def get_paginate_by(self, request):
        rows = request.GET.get('rows', self.get_config(False)['rowNum'])
        try:
            paginate_by = int(rows)
        except ValueError:
            paginate_by = 10
        return paginate_by

    def paginate_df(self, request, df):
        paginate_by = self.get_paginate_by(request)
        if not paginate_by:
            return df, None, None, None

        page_number = int(request.GET.get('page', 1))
        num_records = len(df)
        num_pages = int(math.ceil(num_records / paginate_by))

        start = (page_number - 1) * paginate_by
        end = start + paginate_by + 1
        df_slice = df[start:end]
        return df_slice, page_number, num_pages, num_records

    def get_dataframe_rows(self, df):
        items = []
        for i, row in df.iterrows():
            data = dict(iter(row.items()))
            data['ID'] = i
            items.append(data)
        return items

    def get_json(self, request):
        page_number, num_pages, items, num_records = self.get_items(request)
        data = {
            'page': int(page_number),
            'total': int(num_pages),
            'rows': items,
            'records': int(num_records),
        }
        return json.dumps(data)

    def get_default_config(self):
        config = {
            'datatype': 'json',
            'autowidth': True,
            'forcefit': True,
            'shrinkToFit': True,
            'jsonReader': {'repeatitems': False},
            'rowNum': 10,
            'rowList': [10, 15, 20, 25, 50],
            'sortname': 'ID',
            'viewrecords': True,
            'sortorder': "asc",
            'pager': self.pager_id,
            'altRows': True,
            'gridview': True,
            'height': 'auto',
        }
        return config

    def get_url(self):
        return str(self.url)

    def get_caption(self):
        return self.caption or ''

    def get_config(self, as_json=True):
        config = self.get_default_config()
        config.update(self.extra_config)
        config.update({'url': self.get_url(),
                       'caption': self.get_caption(),
                       'colModel': self.get_colmodels()})
        if as_json:
            config = json.dumps(config)
        return config

    def add_overrides(self, colmodel, column_name):
        override = self.get_override(column_name)
        if override:
            colmodel.update(override)
        return colmodel

    def get_colmodels(self, *args, **kwargs):
        colmodels = []
        df = self.get_dataframe()
        index_colmodel = self.get_index()
        self.add_overrides(index_colmodel, index_colmodel["index"])
        colmodels.append(index_colmodel)
        for column_name in df.columns:
            column = df[column_name]
            colmodel = self.column_to_colmodel(column, column_name)
            self.add_overrides(colmodel, column_name)
            colmodels.append(colmodel)
        return colmodels

    @staticmethod
    def get_index():
        return {
            'name': 'ID',
            'index': 'ID',
            'label': 'ID',
            'editable': False
        }

    @staticmethod
    def column_to_colmodel(column, column_name):
        colmodel = {
            'name': column_name,
            'index': column_name,
            'label': column_name,
            # TODO: dtype?
            'editable': True
        }
        return colmodel
