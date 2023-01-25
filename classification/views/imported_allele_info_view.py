from typing import Optional

from django.contrib.auth.decorators import user_passes_test
from django.db.models import QuerySet, Q, Count, F
from django.http import HttpRequest
from django.shortcuts import render, get_object_or_404
from requests import Response

from classification.models import ImportedAlleleInfo, ImportedAlleleInfoStatus
from genes.hgvs import CHGVS, CHGVSDiff, chgvs_diff_description
from library.guardian_utils import is_superuser
from library.utils import MultiDiff, MultiDiffInput
from snpdb.models import GenomeBuild
from snpdb.views.datatable_view import DatatableConfig, RichColumn, CellData, SortOrder
import re


class ImportedAlleleInfoColumns(DatatableConfig[ImportedAlleleInfo]):

    def render_status(self, data: CellData):
        return ImportedAlleleInfoStatus(data.value).label

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
                key='pk',
                orderable=True,
                default_sort=SortOrder.DESC
            ),
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
                key='latest_validation__include',
                label="Include",
                orderable=True,
                client_renderer='TableFormat.boolean.bind(null, "false_is_error")'
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
        if (transcript_type_not_supported := self.get_query_param('transcript_type_not_supported')) and transcript_type_not_supported == 'true':
            qs = qs.filter(latest_validation__validation_tags__general__transcript_type_not_supported__isnull=False)
        if (transcript_version_change := self.get_query_param('transcript_version_change')) and transcript_version_change == 'true':
            qs = qs.filter(latest_validation__validation_tags__normalize__transcript_version_change__isnull=False) | \
                 qs.filter(latest_validation__validation_tags__liftover__transcript_version_change__isnull=False)
        if (gene_symbol_change := self.get_query_param('gene_symbol_change')) and gene_symbol_change == 'true':
            qs = qs.filter(latest_validation__validation_tags__normalize__gene_symbol_change__isnull=False) | \
                 qs.filter(latest_validation__validation_tags__liftover__gene_symbol_change__isnull=False)
        if (c_nomen_change := self.get_query_param('c_nomen_change')) and c_nomen_change == 'true':
            qs = qs.filter(latest_validation__validation_tags__normalize__c_nomen_change__isnull=False) | \
                 qs.filter(latest_validation__validation_tags__liftover__c_nomen_change__isnull=False)
        if (missing_build := self.get_query_param('missing_build')) and missing_build == 'true':
            # show records where only 1 build is missing
            qs = qs.filter(
                (Q(latest_validation__validation_tags__builds__missing_37__isnull=False) & \
                 Q(latest_validation__validation_tags__builds__missing_38__isnull=True)) |
                (Q(latest_validation__validation_tags__builds__missing_37__isnull=True) & \
                 Q(latest_validation__validation_tags__builds__missing_38__isnull=False))
            )

        if (confirmed := self.get_query_param('confirmed')) and confirmed == 'true':
            qs = qs.filter(latest_validation__confirmed=True)
        if (exclude := self.get_query_param('exclude')) and exclude == 'true':
            qs = qs.filter(latest_validation__include=False)
        if (error_mode := self.get_query_param('errors_mode')) and error_mode == 'true':
            qs = qs.filter(Q(status=ImportedAlleleInfoStatus.FAILED) | Q(grch37__isnull=True) | Q(grch38__isnull=True))

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
    parts = [MultiDiffInput(f"Imported ({allele_info.imported_genome_build_patch_version})", allele_info.imported_c_hgvs)]
    if allele_info.imported_genome_build_patch_version.genome_build == GenomeBuild.grch37():
        parts += [
            MultiDiffInput("Normalised (GRCh37)", allele_info.grch37.c_hgvs if allele_info.grch37 else None),
            MultiDiffInput("Liftover (GRCh38)", allele_info.grch38.c_hgvs if allele_info.grch38 else None)
        ]
    else:
        parts += [
            MultiDiffInput("Normalised (GRCh38)", allele_info.grch38.c_hgvs if allele_info.grch38 else None),
            MultiDiffInput("Liftover (GRCh37)", allele_info.grch37.c_hgvs if allele_info.grch37 else None)
        ]

    diff_output = multi_diff.diffs(parts)

    normalized_diff: Optional[CHGVSDiff] = None
    liftover_diff: Optional[CHGVSDiff] = None
    if imported_c_hgvs := allele_info.imported_c_hgvs_obj:
        if normalized := allele_info.variant_info_for_imported_genome_build:
            if c_hgvs := normalized.c_hgvs_obj:
                normalized_diff = imported_c_hgvs.diff(c_hgvs)
    if (c37 := allele_info.grch37) and (c38 := allele_info.grch38):
        if (c37c := c37.c_hgvs_obj) and (c38c := c38.c_hgvs_obj):
            liftover_diff = c37c.diff(c38c)

    return render(request, "classification/imported_allele_info_detail.html", {
        "allele_info": get_object_or_404(ImportedAlleleInfo, pk=pk),
        "c_hgvses": diff_output,
        "normalized_diff": chgvs_diff_description(normalized_diff) if normalized_diff else None,
        "liftover_diff": chgvs_diff_description(liftover_diff) if liftover_diff else None,
        "variant_coordinate_label": f"Normalised Variant Coordinate ({allele_info.imported_genome_build_patch_version})",
        "validation_tags": allele_info.latest_validation.validation_tags_list if allele_info.latest_validation else None,
        "on_allele_page": request.GET.get("on_allele_page") == "true"
    })
