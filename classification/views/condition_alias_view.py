import urllib

import requests
from django.contrib.auth.models import User
from django.shortcuts import render
from guardian.shortcuts import get_objects_for_user
from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK
import re

from rest_framework.views import APIView

from classification.models import ConditionAlias
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder, BaseDatatableView


ID_EXTRACT_MINI_P = re.compile(r"MONDO:([0-9]+)$")


class ConditionAliasColumns(DatatableConfig):

    def __init__(self, request):
        super().__init__(request)

        self.rich_columns = [
            RichColumn('id', name='ID', client_renderer='renderId', orderable=True),
            RichColumn('lab__name', name='Lab', orderable=True),
            RichColumn('source_gene_symbol', label='Gene Symbol', orderable=True),
            RichColumn('source_text', name='Text', orderable=True),
            RichColumn('records_affected', name='Records Affected', orderable=True, default_sort=SortOrder.DESC, sort_keys=["records_affected", "id"]),
            RichColumn('status', orderable=True),
            RichColumn('created', client_renderer='TableFormat.timestamp', orderable=True)
        ]

    def get_initial_queryset(self):
        return get_objects_for_user(self.user, ConditionAlias.get_read_perm(), klass=ConditionAlias, accept_global_perms=True)


class ConditionAliasDatatableView(BaseDatatableView):

    def config(self, request):
        return ConditionAliasColumns(request)


def condition_aliases_view(request):
    return render(request, 'classification/condition_aliases.html', context={
        'datatable_config': ConditionAliasColumns(request)
    })


def condition_alias_view(request, pk: int):
    user: User = request.user
    condition_alias = ConditionAlias.objects.get(pk=pk)
    condition_alias.check_can_view(user)

    return render(request, 'classification/condition_alias.html', context={
        'condition_alias': condition_alias
    })


class SearchConditionView(APIView):

    def get(self, request, **kwargs) -> Response:
        search_term = request.GET.get('search_term')
        # a regular escape / gets confused for a URL divider
        urllib.parse.quote(search_term).replace('/', '%252F')

        results = requests.get(f'https://api.monarchinitiative.org/api/search/entity/autocomplete/{search_term}', {
            "prefix": "MONDO",
            "rows": 6,
            "minimal_tokenizer": "false",
            "category": "disease"
        }).json().get("docs")

        clean_results = []

        for result in results:
            o_id = result.get('id')
            label = result.get('label')
            if label:
                label = label[len(label) - 1]
            match = result.get('match')
            # if extracted := ID_EXTRACT_MINI_P.match(o_id):
                # id_part = extracted[1]
                # defn = mondo_defns.get(int(id_part))

            highlight = result.get('highlight')

            clean_results.append({
                "id": o_id,
                "label": label,
                "match": match,
                "highlight": highlight
            })

        return Response(status=HTTP_200_OK, data=clean_results)