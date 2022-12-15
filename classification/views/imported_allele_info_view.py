from typing import Optional

from django.contrib.auth.decorators import user_passes_test
from django.db.models import QuerySet, Q, Count
from django.http import HttpRequest
from django.shortcuts import render
from requests import Response

from classification.models import ImportedAlleleInfo
from genes.hgvs import CHGVS
from library.guardian_utils import is_superuser
from snpdb.views.datatable_view import DatatableConfig, RichColumn, CellData


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
        else:
            return {"error": error}


    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.search_box_enabled = True

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
            )
        ]

    def get_initial_queryset(self) -> QuerySet[ImportedAlleleInfo]:
        ae: ImportedAlleleInfo = ImportedAlleleInfo.objects.first()
        return ImportedAlleleInfo.objects.all().annotate(
            classification_count=Count('classification')
        )

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
