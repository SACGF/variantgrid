import urllib
from typing import Optional, Dict, Any, Set

import requests
from django.contrib import messages
from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.shortcuts import render, redirect
from guardian.shortcuts import get_objects_for_user
from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK
import re

from rest_framework.views import APIView

from annotation.models import MonarchDiseaseOntology, MonarchDiseaseOntologyGeneRelationship, \
    MonarchDiseaseOntologyMIMMorbid
from classification.models import ConditionAlias, ConditionAliasJoin, ConditionAliasStatus, Classification
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder, BaseDatatableView


ID_EXTRACT_MINI_P = re.compile(r"MONDO:([0-9]+)$")
# PANEL_APP_OMIM = re.compile(r"MIM#[ ]*([0-9]+)")
PANEL_APP_OMIM = re.compile(r"([0-9]{5,})")

class ConditionAliasColumns(DatatableConfig):

    def pre_render(self, qs: QuerySet):
        all_aliases: Set[int] = set()
        for aliases in qs.values_list('aliases', flat=True):
            if aliases:
                for alias in aliases:
                    try:
                        all_aliases.add(MonarchDiseaseOntology.mondo_id_as_int(alias))
                    except:
                        pass
        mondo: MonarchDiseaseOntology
        for mondo in MonarchDiseaseOntology.objects.filter(pk__in=all_aliases):
            self.mondo_dict[mondo.id_str] = mondo

    def render_aliases(self, row: Dict[str, Any]):
        populated_aliases = []
        if db_aliases := row.get('aliases'):
            for alias in db_aliases:
                pop_alias = {"id": alias}
                populated_aliases.append(pop_alias)

                mondo: MonarchDiseaseOntology
                if mondo := self.mondo_dict.get(alias):
                    pop_alias["name"] = mondo.name

        return {
            'aliases': populated_aliases,
            'join_mode': row.get('join_mode')
        }

    def __init__(self, request):
        super().__init__(request)

        self.mondo_dict: Dict[str, MonarchDiseaseOntology] = dict()

        self.rich_columns = [
            RichColumn('id', name='ID', client_renderer='renderId', orderable=True),
            RichColumn('lab__name', name='Lab', orderable=True),
            RichColumn('source_gene_symbol', label='Gene Symbol', orderable=True),
            RichColumn('source_text', name='Text', orderable=True),
            RichColumn('aliases', name='Aliases', renderer=self.render_aliases, client_renderer='mondo_list'),
            RichColumn('records_affected', name='Records Affected', orderable=True, default_sort=SortOrder.DESC, sort_keys=["records_affected", "id"]),
            RichColumn('status', orderable=True, client_renderer=RichColumn.choices_client_renderer(ConditionAliasStatus.choices))
        ]
        self.extra_columns = ['join_mode']

    def get_initial_queryset(self):
        return get_objects_for_user(self.user, ConditionAlias.get_read_perm(), klass=ConditionAlias, accept_global_perms=True)


class ConditionAliasDatatableView(BaseDatatableView):

    def config_for_request(self, request):
        return ConditionAliasColumns(request)


def condition_aliases_view(request):
    return render(request, 'classification/condition_aliases.html', context={
        'datatable_config': ConditionAliasColumns(request)
    })


def _populateMondoResult(result, gene_symbol) -> Dict:
    if not isinstance(result, dict):
        result = {"id": result}

    mondo_int = MonarchDiseaseOntology.mondo_id_as_int(result.get('id'))
    mondo_record: Optional[MonarchDiseaseOntology]

    url_part: str
    if mondo_record := MonarchDiseaseOntology.objects.filter(pk=mondo_int).first():
        result['definition'] = mondo_record.definition or 'No description provided'
        if 'label' not in result:
            result['label'] = mondo_record.name
        url_part = mondo_record.id_str
        result['id'] = mondo_record.id_str

        relationship: MonarchDiseaseOntologyGeneRelationship
        if relationship := MonarchDiseaseOntologyGeneRelationship.objects.filter(mondo_id=mondo_int, gene_symbol=gene_symbol).first():
            result['relationship'] = relationship.relationship

    else:
        result['definition'] = None

    url_part = result["id"].replace(":", "_")
    result['url'] = f'https://ontology.dev.data.humancellatlas.org/ontologies/mondo/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F{url_part}'

    return result


def _next_condition_alias(user: User, condition_alias: ConditionAlias) -> Optional[ConditionAlias]:
    pending = ConditionAlias.objects.order_by('-records_affected').order_by('id').filter(status=ConditionAliasStatus.PENDING)
    pending = ConditionAlias.filter_for_user(user, pending)
    return pending.first()


def condition_alias_view(request, pk: int):
    user: User = request.user
    condition_alias = ConditionAlias.objects.get(pk=pk)

    if request.method == "GET":
        condition_alias.check_can_view(user)
    else:
        condition_alias.check_can_write(user)
        aliases = request.POST.get('aliases', '')
        if aliases:
            aliases = aliases.split(',')
        else:
            aliases = None

        join_mode = request.POST.get('join_mode')

        condition_alias.aliases = aliases
        condition_alias.join_mode = ConditionAliasJoin(join_mode)
        condition_alias.updated_by = user
        if aliases:
            condition_alias.status = ConditionAliasStatus.RESOLVED
        else:
            if request.POST.get('unmatchable'):
                condition_alias.status = ConditionAliasStatus.UNMATCHABLE
            else:
                condition_alias.status = ConditionAliasStatus.PENDING

        condition_alias.save()

        messages.add_message(request, messages.INFO, message=f"Alias ({pk}) Updated : {condition_alias}")

        next = request.POST.get('next')
        if next:
            if next_ca := _next_condition_alias(user, condition_alias):
                return redirect("condition_alias", pk=next_ca.pk)
            else:
                return redirect("condition_aliases")

        return redirect("condition_alias", pk=pk)

    gene_matches = []
    gene_relationships = MonarchDiseaseOntologyGeneRelationship.objects.filter(gene_symbol=condition_alias.source_gene_symbol)
    mdgr: MonarchDiseaseOntologyGeneRelationship
    for mdgr in gene_relationships:
        gene_match = _populateMondoResult(mdgr.mondo_id, gene_symbol=condition_alias.source_gene_symbol)
        gene_matches.append(gene_match)

    matches_ids = (condition_alias.aliases or [])
    matches = [_populateMondoResult(m_id, gene_symbol=condition_alias.source_gene_symbol) for m_id in matches_ids if m_id]

    panel_app_matches = []
    # try:
    panel_app_matches = search_panelapp(gene_symbol=condition_alias.source_gene_symbol)
    # except:
    #    pass

    return render(request, 'classification/condition_alias.html', context={
        'classification_subset': condition_alias.classification_modifications[0:4],
        'extra_count': max(condition_alias.classification_modifications.count() - 4, 0),
        'condition_alias': condition_alias,
        'matches': matches,
        'panel_app_matches': panel_app_matches,
        'gene_matches': gene_matches
    })


def search_panelapp(gene_symbol: str):
    results = requests.get(f'https://panelapp.agha.umccr.org/api/v1/genes/{gene_symbol}/').json().get("results")
    mondo_results = {}
    for panel_app_result in results:
        phenotype_row: str
        for phenotype_row in panel_app_result.get("phenotypes", []):
            for omim_match in PANEL_APP_OMIM.finditer(phenotype_row):
                omim_id = omim_match[1]
                if mondo_rel := MonarchDiseaseOntologyMIMMorbid.objects.filter(omim_id=omim_id).first():
                    mondo = mondo_rel.mondo
                    if not mondo in mondo_results:
                        mondo_result = _populateMondoResult(mondo.pk, gene_symbol)
                        mondo_results[mondo.pk] = mondo_result

                        mondo_result["panelapp_evidence"] = set()
                        mondo_result["panelapp_phenotype"] = phenotype_row
                        evidence = mondo_result.get("panelapp_evidence")
                    for evidence_str in panel_app_result.get("evidence"):
                        evidence.add(evidence_str)

    for mondo_result in mondo_results.values():
        mondo_result["panelapp_evidence"] = list(mondo_result["panelapp_evidence"])
        mondo_result["panelapp_evidence"].sort()

    return list(mondo_results.values())

class SearchConditionView(APIView):

    def get(self, request, **kwargs) -> Response:
        search_term = request.GET.get('search_term')
        gene_symbol = request.GET.get('gene_symbol')
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
                label = label[0]
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

        # populate more with our database
        # do this separate so if we cache the above we can still apply the below
        for result in clean_results:
            _populateMondoResult(result, gene_symbol=gene_symbol)

        return Response(status=HTTP_200_OK, data=clean_results)
