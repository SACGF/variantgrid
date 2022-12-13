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
        build_number = '37'
        if '38' in data.key:
            build_number = '38'
        if c_hgvs_str := data.get(f'grch{build_number}__c_hgvs'):
            c_hgvs = CHGVS(c_hgvs_str)
            variant_id = data.get(f'grch{build_number}__variant')

            return {
                "transcript": c_hgvs.transcript,
                "geneSymbol": c_hgvs.gene_symbol,
                "variant": c_hgvs.variant,
                "variantId": variant_id
            }
        else:
            return {}


    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.search_box_enabled = True

        self.rich_columns = [
            RichColumn(
                key='imported_c_hgvs',
                label='Imported c.HGVS',
                orderable=True
            ),
            RichColumn(
                key='imported_genome_build',
                label='Imported Gene Symbol',
                orderable=True
            ),
            RichColumn(
                key='grch37__c_hgvs',
                label='Resolved GRCh37 c.HGVS',
                orderable=True,
                extra_columns=['grch37__variant'],
                renderer=self.render_c_hgvs,
                client_renderer='TableFormat.hgvs'
            ),
            RichColumn(
                key='grch38__c_hgvs',
                label='Resolved GRCh38 c.HGVS',
                orderable=True,
                extra_columns=['grch38__variant'],
                renderer=self.render_c_hgvs,
                client_renderer='TableFormat.hgvs'
            ),
            RichColumn(
                key='classification_count',
                label="Classification Count",
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
