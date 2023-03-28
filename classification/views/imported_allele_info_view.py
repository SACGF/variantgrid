import operator
import re
from functools import reduce
from typing import Optional

from django.db.models import QuerySet, Q, Count
from django.http import HttpRequest
from django.shortcuts import render, get_object_or_404
from requests import Response

from classification.models import ImportedAlleleInfo, ImportedAlleleInfoStatus, Classification
from classification.models.classification_variant_info_models import ImportedAlleleInfoValidation
from genes.hgvs import CHGVS, CHGVSDiff, chgvs_diff_description
from library.django_utils import get_url_from_view_path
from library.utils import MultiDiff, MultiDiffInput, ExportRow, export_column
from library.utils.django_utils import render_ajax_view
from snpdb.admin_utils import get_admin_url
from snpdb.models import GenomeBuild, Lab, Allele
from snpdb.views.datatable_view import DatatableConfig, RichColumn, CellData, SortOrder


class ImportedAlleleInfoColumns(DatatableConfig[ImportedAlleleInfo]):

    @staticmethod
    def render_c_hgvs(data: CellData):
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
                # json_data['variant_id'] = variant_id
                return json_data
            else:
                return {
                    "full": c_hgvs_str
                }
        elif data.key == 'imported_c_hgvs':
            if g_hgvs := data.get('imported_g_hgvs'):
                return {"full": f'(g.HGVS) {g_hgvs}'}
            else:
                return {"full": None}
        elif not error:
            if not data.get('variant_coordinate'):
                return {"error": "Could not derive variant coordinates"}
            elif variant_id:
                return {"error": "Resolved to variant, but no c.hgvs generated"}
            else:
                return {"error": "Not resolved to a variant"}
        else:
            return {"error": error}

    @staticmethod
    def render_allele(data: CellData):
        if value := data.value:
            allele = Allele.objects.get(pk=value)
            return {
                "text": f'{allele:CA}',
                "url": allele.get_absolute_url()
            }

    @staticmethod
    def render_validation(data: CellData):
        return {
            "include": data.get('latest_validation__include') is True,
            "tags": [tag.as_json() for tag in ImportedAlleleInfoValidation.validation_tags_list_from_dict(data.get('latest_validation__validation_tags'))]
        }

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
                renderer=ImportedAlleleInfoColumns.render_c_hgvs,
                client_renderer='VCTable.hgvs',
                extra_columns=['imported_g_hgvs', 'variant_coordinate']
            ),
            RichColumn(
                key='grch37__c_hgvs',
                label='Resolved GRCh37<br/>c.HGVS',
                orderable=True,
                extra_columns=['grch37__variant', 'grch37__error'],
                renderer=ImportedAlleleInfoColumns.render_c_hgvs,
                client_renderer='VCTable.hgvs'
            ),
            RichColumn(
                key='grch38__c_hgvs',
                label='Resolved GRCh38<br/>c.HGVS',
                orderable=True,
                extra_columns=['grch38__variant', 'grch38__error'],
                renderer=ImportedAlleleInfoColumns.render_c_hgvs,
                client_renderer='VCTable.hgvs'
            ),
            # RichColumn(
            #     key='allele',
            #     label='Allele',
            #     orderable=True,
            #     renderer=ImportedAlleleInfoColumns.render_allele,
            #     client_renderer='TableFormat.linkUrl'
            # ),
            RichColumn(
                key='latest_validation__include',
                label="Include",
                orderable=True,
                extra_columns=['latest_validation__validation_tags'],
                renderer=ImportedAlleleInfoColumns.render_validation,
                client_renderer='render_validation'
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
                (Q(latest_validation__validation_tags__builds__missing_37__isnull=False) &
                 Q(latest_validation__validation_tags__builds__missing_38__isnull=True)) |
                (Q(latest_validation__validation_tags__builds__missing_37__isnull=True) &
                 Q(latest_validation__validation_tags__builds__missing_38__isnull=False))
            )

        if (confirmed := self.get_query_param('confirmed')) and confirmed == 'true':
            qs = qs.filter(latest_validation__confirmed=True)
        if (exclude := self.get_query_param('exclude')) and exclude == 'true':
            qs = qs.filter(latest_validation__include=False)
        if (error_mode := self.get_query_param('errors_mode')) and error_mode == 'true':
            qs = qs.filter(Q(status=ImportedAlleleInfoStatus.FAILED) | Q(latest_validation__validation_tags__general__hgvs_issue__isnull=False))  # | Q(grch37__isnull=True) | Q(grch38__isnull=True))
        if (in_progress := self.get_query_param('in_progress')) and in_progress == 'true':
            qs = qs.filter(Q(status__in=(ImportedAlleleInfoStatus.PROCESSING, ImportedAlleleInfoStatus.MATCHED_IMPORTED_BUILD)))

        return qs

    def power_search(self, qs: QuerySet[ImportedAlleleInfo], search_string: str) -> QuerySet[ImportedAlleleInfo]:
        # TODO, make RichColumn's searchable on/off so we can just fall back onto that
        ors = [
            Q(imported_c_hgvs__icontains=search_string),
            Q(grch37__c_hgvs=search_string),
            Q(grch38__c_hgvs=search_string)
        ]
        try:
            id_int = int(search_string)
            ors.append(Q(pk=id_int))
        except ValueError:
            pass
        return qs.filter(reduce(operator.or_, ors))


def view_imported_allele_info(request: HttpRequest) -> Response:
    return render(request, "classification/imported_allele_info.html", {})

def view_imported_allele_info_detail(request: HttpRequest, allele_info_id: int):
    allele_info = get_object_or_404(ImportedAlleleInfo, pk=allele_info_id)
    # just split up c.hgvs into logical parts, and then the diff will reset with each new group (treat it as different words)

    HGVS_BASE_REGEX_STR = '^(?P<transcript>[^.]+?)(?P<transcript_version>\.[0-9]+)?(?P<gene_symbol>[(].*[)])?(?P<c_dot>:[cng]\.)'

    HGVS_REGEX_BASIC = re.compile(HGVS_BASE_REGEX_STR + '(?P<c_nomen_pos>[0-9+_-]*)(?P<c_nomen_change>.*?)$')
    HGVS_REGEX_REF_ALT = re.compile(HGVS_BASE_REGEX_STR + '(?P<c_nomen_pos>[0-9+_-]*?)(?P<ref>[ACTG]+)(?P<operation>>)(?P<alt>[ACTG]+)$')
    HGVS_REGEX_DEL_INS = re.compile(HGVS_BASE_REGEX_STR + '(?P<c_nomen_pos>[0-9+_-]*?)(?P<del>del)(?P<ref>[ACTG]*)(?P<ins>del)(?P<alt>[ACTG]*)$')
    HGVS_REGEX_SIMPLE_OP = re.compile(HGVS_BASE_REGEX_STR + '(?P<c_nomen_pos>[0-9+_-]*?)(?P<operation>dup|del|ins)(?P<alt>[ACTG]*)$')

    use_regex = HGVS_REGEX_BASIC
    regex_attempt_order = [HGVS_REGEX_REF_ALT, HGVS_REGEX_DEL_INS, HGVS_REGEX_SIMPLE_OP]
    for regex in regex_attempt_order:
        if regex.match(allele_info.imported_c_hgvs):
            use_regex = regex
            break

    multi_diff = MultiDiff(use_regex)
    parts = [MultiDiffInput(f"Imported ({allele_info.imported_genome_build_patch_version})", allele_info.imported_c_hgvs)]
    if allele_info.imported_genome_build_patch_version.genome_build == GenomeBuild.grch37():
        parts += [
            MultiDiffInput("Normalised (GRCh37)", allele_info.grch37.c_hgvs if allele_info.grch37 else None, is_reference=True),
            MultiDiffInput("Liftover (GRCh38)", allele_info.grch38.c_hgvs if allele_info.grch38 else None)
        ]
    else:
        parts += [
            MultiDiffInput("Normalised (GRCh38)", allele_info.grch38.c_hgvs if allele_info.grch38 else None, is_reference=True),
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

    return render_ajax_view(request, "classification/imported_allele_info_detail.html", {
        "allele_info": allele_info,
        "c_hgvses": diff_output,
        "normalized_diff": chgvs_diff_description(normalized_diff) if normalized_diff else None,
        "liftover_diff": chgvs_diff_description(liftover_diff) if liftover_diff else None,
        "variant_coordinate_label": f"Normalised Variant Coordinate ({allele_info.imported_genome_build_patch_version})",
        "validation_tags": allele_info.latest_validation.validation_tags_list if allele_info.latest_validation else None,
        "on_allele_page": request.GET.get("on_allele_page") == "true"
    }, menubar='classification')


class ImportedAlleleInfoDownload(ExportRow):

    def __init__(self, allele_info: ImportedAlleleInfo):
        self.allele_info = allele_info

    @export_column(label="ID")
    def id(self):
        return self.allele_info.id

    @export_column("URL")
    def url(self):
        return get_url_from_view_path(get_admin_url(self.allele_info))

    @export_column(label="Imported Genome Build")
    def imported_genome_build(self):
        return str(self.allele_info.imported_genome_build_patch_version)

    @export_column(label="c.HGVS (Imported)")
    def c_hgvs_imported(self):
        return self.allele_info.imported_c_hgvs

    @export_column(label="c.HGVS (37)")
    def c_hgvs_37(self):
        if c37 := self.allele_info[GenomeBuild.grch37()]:
            return c37.c_hgvs

    @export_column(label="c.HGVS (38)")
    def c_hgvs_38(self):
        if c38 := self.allele_info[GenomeBuild.grch37()]:
            return c38.c_hgvs

    @export_column(label="Differences")
    def differences(self):
        return "\n".join(str(tag) for tag in self.allele_info.latest_validation.validation_tags_list)

    @export_column(label="Included")
    def included(self):
        return self.allele_info.latest_validation.include

    @export_column(label="Confirmed")
    def confirmed(self):
        return self.allele_info.latest_validation.confirmed

    @export_column(label="Classification Count")
    def classification_count(self):
        return Classification.objects.filter(allele_info=self.allele_info).count()

    @export_column(label="Allele URL")
    def classification_count(self):
        if allele := self.allele_info.allele:
            return get_url_from_view_path(allele.get_absolute_url())

    @export_column(label="Involved Labs")
    def involved_labs(self):
        return ", ".join([str(lab) for lab in sorted(Lab.objects.filter(pk__in=Classification.objects.filter(allele_info=self.allele_info).values_list('lab', flat=True)).select_related('organization'))])


def download_allele_info(request: HttpRequest):
    return ImportedAlleleInfoDownload.streaming_csv(
        ImportedAlleleInfo.objects.order_by('id').iterator(chunk_size=2000),
        filename="imported_allele_infos"
    )
