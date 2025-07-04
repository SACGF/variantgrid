import operator
import re
from functools import reduce
from typing import Optional

from django.db.models import QuerySet, Q, Count
from django.http import HttpRequest
from django.shortcuts import render, get_object_or_404
from requests import Response

from classification.models import ImportedAlleleInfo, ImportedAlleleInfoStatus, Classification, \
    ClassificationModification
from classification.models.classification_variant_info_models import ImportedAlleleInfoValidation
from genes.hgvs import CHGVS, CHGVSDiff, chgvs_diff_description
from library.django_utils import get_url_from_view_path
from library.django_utils import require_superuser
from library.utils import MultiDiff, MultiDiffInput, ExportRow, export_column
from library.utils.django_utils import render_ajax_view
from snpdb.admin_utils import get_admin_url
from snpdb.models import GenomeBuild, Lab, Allele
from snpdb.views.datatable_view import DatatableConfig, RichColumn, CellData, SortOrder


class ImportedAlleleInfoColumns(DatatableConfig[ImportedAlleleInfo]):

    @staticmethod
    def render_c_hgvs(data: CellData):
        c_hgvs_str: Optional[str]
        c_hgvs: Optional[CHGVS]
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
            "tags": [tag.as_json() for tag in ImportedAlleleInfoValidation.validation_tags_list_from_dict(data.get('latest_validation__validation_tags'))],
            "message": data.get('message')
        }

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.search_box_enabled = True
        self.expand_client_renderer = DatatableConfig._row_expand_ajax('view_imported_allele_info_detail', expected_height=108)

        self.rich_columns = [
            RichColumn(
                key='pk',
                orderable=True,
                default_sort=SortOrder.DESC,
                order_sequence=[SortOrder.DESC, SortOrder.ASC]
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
            RichColumn(
                key='latest_validation__include',
                label="Include",
                orderable=True,
                extra_columns=['latest_validation__validation_tags', 'message'],
                renderer=ImportedAlleleInfoColumns.render_validation,
                client_renderer='render_validation'
            ),
            RichColumn(
                key='classification_count',
                label="<span style=\"font-size:10px\">Classification<br/>Record Count</span>",
                order_sequence=[SortOrder.DESC, SortOrder.ASC]
            ),
            RichColumn('id', visible=False)  # just used for the expand
        ]

    def get_initial_queryset(self) -> QuerySet[ImportedAlleleInfo]:
        return ImportedAlleleInfo.objects.all().annotate(
            classification_count=Count('classification', filter=Q(classification__withdrawn=False))
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

        if self.get_query_param('confirmed') == 'true':
            qs = qs.filter(latest_validation__confirmed=True)
        if self.get_query_param('exclude') == 'true':
            qs = qs.filter(latest_validation__include=False)
        if self.get_query_param('errors_mode') == 'true':
            qs = qs.filter(Q(status=ImportedAlleleInfoStatus.FAILED) | Q(latest_validation__validation_tags__general__hgvs_issue__isnull=False))  # | Q(grch37__isnull=True) | Q(grch38__isnull=True))
        if self.get_query_param('in_progress') == 'true':
            qs = qs.filter(Q(status__in=(ImportedAlleleInfoStatus.PROCESSING, ImportedAlleleInfoStatus.MATCHED_IMPORTED_BUILD)))
        if self.get_query_param('dirty') == 'true':
            qs = qs.filter(Q(dirty_message__isnull=False))

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


@require_superuser
def view_imported_allele_info(request: HttpRequest) -> Response:
    status = request.GET.get("status")
    return render(request, "classification/imported_allele_info.html", {"status": status})


def view_imported_allele_info_detail(request: HttpRequest, allele_info_id: int):
    allele_info = get_object_or_404(ImportedAlleleInfo, pk=allele_info_id)
    # just split up c.hgvs into logical parts, and then the diff will reset with each new group (treat it as different words)

    HGVS_BASE_REGEX_STR = r'^(?P<transcript>[^.]+?)(?P<transcript_version>\.[0-9]+)?(?P<gene_symbol>[(].*[)])?(?P<c_dot>:[cng]\.)'

    HGVS_REGEX_BASIC = re.compile(HGVS_BASE_REGEX_STR + '(?P<c_nomen_full>.*?)$')
    HGVS_REGEX_REF_ALT = re.compile(HGVS_BASE_REGEX_STR + '(?P<c_nomen_pos>[0-9+_-]*?)(?P<ref>[ACTG]+)(?P<operation>>)(?P<alt>[ACTG]+)$')
    HGVS_REGEX_DEL_INS = re.compile(HGVS_BASE_REGEX_STR + '(?P<c_nomen_pos>[0-9+_-]*?)(?P<del>del)(?P<ref>[ACTG]*)(?P<ins>ins)(?P<alt>[ACTG]*)$')
    HGVS_REGEX_SIMPLE_OP = re.compile(HGVS_BASE_REGEX_STR + '(?P<c_nomen_pos>[0-9+_-]*?)(?P<operation>dup|del|ins)(?P<alt>[ACTG]*)$')

    FALLBACK_HGVS = re.compile("(?P<all>.*)")
    use_texts = [text for text in [allele_info.imported_c_hgvs, allele_info.grch37.c_hgvs if allele_info.grch37 else None, allele_info.grch38.c_hgvs if allele_info.grch38 else None] if text]

    use_regex = FALLBACK_HGVS
    if use_texts:
        regex_attempt_order = [HGVS_REGEX_REF_ALT, HGVS_REGEX_DEL_INS, HGVS_REGEX_SIMPLE_OP, HGVS_REGEX_BASIC]
        for regex in regex_attempt_order:
            is_a_match = True
            for use_text in use_texts:
                if not regex.match(use_text):
                    is_a_match = False
                    break

            if is_a_match:
                use_regex = regex
                break

    multi_diff = MultiDiff(use_regex)
    # 'parts' and 'resolved_variant_info' are kept in sync
    parts = []
    c_hgvs_resolved_variant_info = []

    if allele_info.imported_c_hgvs:
        parts += [MultiDiffInput(f"Imported ({allele_info.imported_genome_build_patch_version})", allele_info.imported_c_hgvs)]
        c_hgvs_resolved_variant_info.append(None)

    if allele_info.imported_genome_build_patch_version.genome_build == GenomeBuild.grch37():
        resolved_variant_info_is_reference = [
            (allele_info.grch37, True),
            (allele_info.grch38, False),
        ]
    else:
        resolved_variant_info_is_reference = [
            (allele_info.grch38, True),
            (allele_info.grch37, False),
        ]

    for rvi, is_reference in resolved_variant_info_is_reference:
        if is_reference:
            origin = "Normalized"
        else:
            origin = "Liftover"
        if rvi:
            label = f"{origin} ({rvi.genome_build})"
        else:
            label = origin

        parts.append(MultiDiffInput(label, rvi.c_hgvs if rvi else None,
                                    is_reference=is_reference))
        c_hgvs_resolved_variant_info.append((label + " c.HGVS by", rvi.c_hgvs_converter_version if rvi else ""))

    diff_output = multi_diff.diffs(parts)
    if not allele_info.imported_c_hgvs:
        diff_output = [None] + diff_output

    normalized_diff: Optional[CHGVSDiff] = None
    liftover_diff: Optional[CHGVSDiff] = None
    if imported_c_hgvs := allele_info.imported_c_hgvs_obj:
        if normalized := allele_info.variant_info_for_imported_genome_build:
            if c_hgvs := normalized.c_hgvs_obj:
                normalized_diff = imported_c_hgvs.diff(c_hgvs)
    if (c37 := allele_info.grch37) and (c38 := allele_info.grch38):
        if (c37c := c37.c_hgvs_obj) and (c38c := c38.c_hgvs_obj):
            liftover_diff = c37c.diff(c38c)

    classifications = ClassificationModification.latest_for_user(user=request.user, published=True).filter(classification__in=allele_info.classification_set.all())

    return render_ajax_view(request, "classification/imported_allele_info_detail.html", {
        "allele_info": allele_info,
        "g_hgvs_label": f"Imported g.HGVS ({allele_info.imported_genome_build_patch_version})",
        "c_hgvses": diff_output,
        "c_hgvs_resolved_variant_info": c_hgvs_resolved_variant_info,
        "normalized_diff": chgvs_diff_description(normalized_diff) if normalized_diff else None,
        "liftover_diff": chgvs_diff_description(liftover_diff) if liftover_diff else None,
        "variant_coordinate_label": f"Variant Coordinate (from HGVS resolution) ({allele_info.imported_genome_build_patch_version})",
        "variant_coordinate_normalized_label": f"Variant Coordinate Normalized ({allele_info.imported_genome_build_patch_version})",
        "validation_tags": allele_info.latest_validation.validation_tags_list if allele_info.latest_validation else None,
        "on_allele_page": request.GET.get("on_allele_page") == "true",
        "classifications": classifications
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
        """"
        Technically should be called "HGVS (Imported)" but leaving as c.HGVS for compatibility with tools
        """
        if imported_c_hgvs := self.allele_info.imported_c_hgvs:
            return imported_c_hgvs
        else:
            return self.allele_info.imported_g_hgvs

    @export_column(label="c.HGVS (37)")
    def c_hgvs_37(self):
        if c37 := self.allele_info[GenomeBuild.grch37()]:
            return c37.c_hgvs

    @export_column(label="c.HGVS (38)")
    def c_hgvs_38(self):
        if c38 := self.allele_info[GenomeBuild.grch38()]:
            return c38.c_hgvs

    @export_column(label="Differences")
    def differences(self):
        return "\n".join(str(tag) for tag in self.allele_info.latest_validation.validation_tags_list)

    @export_column(label="c.HGVS (37) sort")
    def c_hgvs_37_sort(self):
        if c37 := self.allele_info[GenomeBuild.grch37()]:
            return c37.genomic_sort

    @export_column(label="c.HGVS (38) sort")
    def c_hgvs_38_sort(self):
        if c38 := self.allele_info[GenomeBuild.grch38()]:
            return c38.genomic_sort

    @export_column(label="variant ID (37)")
    def variant_id_37(self):
        if resolved := self.allele_info[GenomeBuild.grch37()]:
            return resolved.variant_id

    @export_column(label="variant ID (38)")
    def variant_id_38(self):
        if resolved := self.allele_info[GenomeBuild.grch38()]:
            return resolved.variant_id

    @export_column(label="Included")
    def included(self):
        return self.allele_info.latest_validation.include

    @export_column(label="Confirmed")
    def confirmed(self):
        return self.allele_info.latest_validation.confirmed

    @export_column(label="Classification Record Count")
    def classification_count(self):
        return Classification.objects.filter(allele_info=self.allele_info, withdrawn=False).count()

    @export_column(label="Allele URL")
    def allele_url(self):
        if allele := self.allele_info.allele:
            return get_url_from_view_path(allele.get_absolute_url())

    @export_column(label="Involved Labs")
    def involved_labs(self):
        return ", ".join([str(lab) for lab in sorted(Lab.objects.filter(pk__in=Classification.objects.filter(allele_info=self.allele_info, withdrawn=False).values_list('lab', flat=True)).select_related('organization'))])


@require_superuser
def download_allele_info(request: HttpRequest):
    return ImportedAlleleInfoDownload.streaming_csv(
        ImportedAlleleInfo.objects.order_by('id').iterator(chunk_size=2000),
        filename="imported_allele_infos"
    )
