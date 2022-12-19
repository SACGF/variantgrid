import itertools
from dataclasses import dataclass, field
from enum import Enum
from itertools import zip_longest
from typing import Optional, List, Pattern, Any

from django.contrib.auth.decorators import user_passes_test
from django.db.models import QuerySet, Q, Count, F
from django.http import HttpRequest
from django.shortcuts import render, get_object_or_404
from requests import Response

from classification.models import ImportedAlleleInfo
from genes.hgvs import CHGVS
from library.guardian_utils import is_superuser
from library.utils import first, MultiDiff, MultiDiffInput
from snpdb.views.datatable_view import DatatableConfig, RichColumn, CellData
import re


class ImportedAlleleInfoColumns(DatatableConfig[ImportedAlleleInfo]):

    def render_c_hgvs(self, data: CellData):
        c_hgvs_str: Optional[str] = None
        c_hgvs: Optional[CHGVS] = None
        variant_id: Optional[int] = None
        error: Optional[str] = None

        if '37' in data.key:
            c_hgvs_str = data.get('grch37__c_hgvs')
            variant_id = data.get('grch37__variant')
            error = data.get('grch37__error')
        elif '38' in data.key:
            c_hgvs_str = data.get('grch38__c_hgvs')
            variant_id = data.get('grch38__variant')
            error = data.get('grch38__error')
        else:
            c_hgvs_str = data.get(data.key)

        if c_hgvs_str:
            if c_hgvs := CHGVS(c_hgvs_str):
                json_data = c_hgvs.to_json()
                json_data['variant_id'] = variant_id
                return json_data
            else:
                return {
                    "full": c_hgvs_str
                }
        elif not variant_id:
            return {"error": "Not resolved to a variant"}
        else:
            return {"error": error}


    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.search_box_enabled = True
        self.expand_client_renderer = DatatableConfig._row_expand_ajax('view_imported_allele_info_detail', expected_height=108)

        self.rich_columns = [
            RichColumn(
                key='imported_genome_build_patch_version',
                label='Imported<br/>Genome Build',
                orderable=True
            ),
            RichColumn(
                key='imported_c_hgvs',
                label='Imported<br/>c.HGVS',
                orderable=True,
                renderer=self.render_c_hgvs,
                client_renderer='VCTable.hgvs'
            ),
            RichColumn(
                key='grch37__c_hgvs',
                label='Resolved GRCh37<br/>c.HGVS',
                orderable=True,
                extra_columns=['grch37__variant', 'grch37__error'],
                renderer=self.render_c_hgvs,
                client_renderer='VCTable.hgvs'
            ),
            RichColumn(
                key='grch38__c_hgvs',
                label='Resolved GRCh38<br/>c.HGVS',
                orderable=True,
                extra_columns=['grch38__variant', 'grch38__error'],
                renderer=self.render_c_hgvs,
                client_renderer='VCTable.hgvs'
            ),
            RichColumn(
                key='classification_count',
                label="Classification<br/>Count",
                orderable=True
            ),
            RichColumn('id', visible=False)  # just used for the expand
        ]

    def get_initial_queryset(self) -> QuerySet[ImportedAlleleInfo]:
        ae: ImportedAlleleInfo = ImportedAlleleInfo.objects.first()
        return ImportedAlleleInfo.objects.all().annotate(
            classification_count=Count('classification')
        )

    def filter_queryset(self, qs: QuerySet[ImportedAlleleInfo]) -> QuerySet[ImportedAlleleInfo]:
        if filter_mode := self.get_query_param('37_not_38'):
            if filter_mode == 'true':
                qs = qs.filter(grch37__c_hgvs__ne=F('grch38__c_hgvs'))
        return qs

    def power_search(self, qs: QuerySet[ImportedAlleleInfo], search_string: str) -> QuerySet[ImportedAlleleInfo]:
        # TODO, make RichColumn's searchable on/off so we can just fall back onto that
        return qs.filter(
            Q(imported_c_hgvs__icontains=search_string) | \
            Q(grch37__c_hgvs=search_string) | \
            Q(grch38__c_hgvs=search_string)
        )

@user_passes_test(is_superuser)
def view_imported_allele_info(request: HttpRequest) -> Response:
    return render(request, "classification/imported_allele_info.html", {})


@user_passes_test(is_superuser)
def view_imported_allele_info_detail(request: HttpRequest, pk: int):
    allele_info = get_object_or_404(ImportedAlleleInfo, pk=pk)
    # just split up c.hgvs into logical parts, and then the diff will reset with each new group (treat it as different words)
    HGVS_REGEX = re.compile(
        '(?P<transcript>[^.]+?)'
        '(?P<transcript_version>\.[0-9]+)?'
        '(?P<gene_symbol>[(].*[)])?'
        '(?P<c_dot>:c\.)'
        '(?P<c_nomen_pos>[0-9]+)'
        '(?P<c_nomen_change>.*)'
    )
    multi_diff = MultiDiff(HGVS_REGEX)
    parts = [
        MultiDiffInput(f"Imported ({allele_info.imported_genome_build_patch_version})", allele_info.imported_c_hgvs),
        MultiDiffInput("GRCh37", allele_info.grch37.c_hgvs if allele_info.grch37 else None),
        MultiDiffInput("GRCh38", allele_info.grch38.c_hgvs if allele_info.grch38 else None)
    ]
    diff_output = multi_diff.diffs(parts)

    return render(request, "classification/imported_allele_info_detail.html", {
        "allele_info": get_object_or_404(ImportedAlleleInfo, pk=pk),
        "c_hgvses": diff_output
    })