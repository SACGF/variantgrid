import json
from collections import defaultdict
from typing import Dict, Any, Optional, Iterable

from django.db.models import QuerySet, When, Value, Case, IntegerField, Count, Q
from django.http import HttpResponse, StreamingHttpResponse, HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render, redirect, get_object_or_404
from django.urls import reverse
from django.views import View
from lazy import lazy

from classification.enums import SpecialEKeys
from classification.models import ClinVarExport, ClinVarExportBatch, ClinVarExportBatchStatus, \
    EvidenceKeyMap, ClinVarExportStatus, ClinVarExportSubmission
from classification.utils.clinvar_matcher import ClinVarLegacyRow
from classification.views.classification_dashboard_view import ClassificationDashboard
from genes.hgvs import CHGVS
from library.cache import timed_cache
from library.django_utils import add_save_message, get_url_from_view_path
from library.utils import html_to_text, export_column, ExportRow, local_date_string
from snpdb.models import ClinVarKey, Lab, Allele, GenomeBuild
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder
from uicore.json.json_types import JsonDataType
import io


@timed_cache(size_limit=30, ttl=60)
def allele_for(allele_id: int) -> Allele:
    return Allele.objects.select_related('clingen_allele').get(pk=allele_id)


class ClinVarExportBatchColumns(DatatableConfig):

    def render_status(self, row: Dict[str, Any]):
        return ClinVarExportBatchStatus(row['status']).label

    def __init__(self, request):
        super().__init__(request)

        self.expand_client_renderer = DatatableConfig._row_expand_ajax('clinvar_export_batch_detail', expected_height=120)
        self.rich_columns = [
            RichColumn("id", label="ID", orderable=True, default_sort=SortOrder.DESC),
            RichColumn("clinvar_key", label="ClinVar Key", orderable=True, enabled=False),
            RichColumn("created", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn("record_count", label="Records In Submission", orderable=True),
            RichColumn("status", renderer=self.render_status, orderable=True)
        ]

    def get_initial_query_params(self, clinvar_key: Optional[str] = None):
        clinvar_keys = ClinVarKey.clinvar_keys_for_user(self.user)
        if clinvar_key:
            clinvar_keys = clinvar_keys.filter(pk=clinvar_key)
        cve = ClinVarExportBatch.objects.filter(clinvar_key__in=clinvar_keys)
        return cve

    def get_initial_queryset(self):
        initial_qs = self.get_initial_query_params(
            clinvar_key=self.get_query_param('clinvar_key')
        )
        initial_qs = initial_qs.annotate(record_count=Count('clinvarexportsubmission'))
        return initial_qs


def clinvar_export_batch_detail(request, clinvar_export_batch_id: int):
    clinvar_export_batch: ClinVarExportBatch = ClinVarExportBatch.objects.get(pk=clinvar_export_batch_id)
    clinvar_export_batch.clinvar_key.check_user_can_access(request.user)
    submissions = clinvar_export_batch.clinvarexportsubmission_set.order_by('created')

    return render(request, 'classification/clinvar_export_batch_detail.html', {
        "batch": clinvar_export_batch,
        "submissions": submissions
    })


class ClinVarExportColumns(DatatableConfig[ClinVarExport]):

    def render_allele(self, row: Dict[str, Any]) -> str:
        allele = allele_for(row["clinvar_allele__allele"])
        return f"{allele:CA}"

    def pre_render(self, qs: QuerySet[ClinVarExport]):
        # find all the batches these records are in
        # do this once rather than per row
        super().pre_render(qs)
        export_to_batches = defaultdict(list)
        for export_id, batch_id in ClinVarExportSubmission.objects.filter(clinvar_export__in=qs.values_list("id", flat=True)).order_by(
            "-submission_batch").values_list("clinvar_export", "submission_batch"):
            export_to_batches[export_id].append(batch_id)
        self.export_to_batches = export_to_batches

    def power_search(self, qs: QuerySet[ClinVarExport], search_string: str) -> QuerySet[ClinVarExport]:
        try:
            batch_id = int(search_string)
            submissions = ClinVarExportSubmission.objects.filter(submission_batch_id=batch_id).values_list("clinvar_export_id", flat=True)
            # filter rather than just return submissions to make sure user has permission to see items that belong to the batch
            # also make sure an individual export doesn't have the batch ID
            return qs.filter(Q(pk__in=submissions) | Q(pk=batch_id))
        except ValueError:
            pass
        return super().power_search(qs, search_string)

    def batches(self, row: Dict[str, Any]) -> JsonDataType:
        return self.export_to_batches.get(row.get("id"))

    def render_c_hgvs(self, row: Dict[str, Any]) -> JsonDataType:
        if row["classification_based_on__classification__variant"]:
            genome_build = row["classification_based_on__published_evidence__genome_build__value"]
            c_hgvs_str: str
            if "h37" in genome_build:
                c_hgvs_str = row["classification_based_on__classification__chgvs_grch37"]
            else:
                c_hgvs_str = row["classification_based_on__classification__chgvs_grch38"]

            data: Dict[str, Any]
            c_hgvs = CHGVS(c_hgvs_str)
            if c_hgvs.raw_c != c_hgvs.full_c_hgvs:
                data = {
                    "genomeBuild": genome_build,
                    "transcript": c_hgvs.transcript,
                    "geneSymbol": c_hgvs.gene_symbol,
                    "variant": c_hgvs.raw_c,
                    # below row makes this link to the variant, but probably not desired action
                    # "variantId": row["classification_based_on__classification__variant"]
                }
            else:
                data = {"full": c_hgvs_str}

            if allele_id := row["clinvar_allele__allele"]:
                allele = allele_for(allele_id)
                data["allele"] = f"{allele:CA}"

            return data
        else:
            return None

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.export_to_batches = dict()

        self.search_box_enabled = True
        self.expand_client_renderer = DatatableConfig._row_expand_ajax('clinvar_export_detail')
        self.rich_columns = [
            RichColumn("id", label="ID", orderable=True, default_sort=SortOrder.DESC),
            RichColumn(name="c_hgvs", label='Allele',
                    sort_keys=["classification_based_on__classification__chgvs_grch38"],
                    extra_columns=[
                        "clinvar_allele__allele",
                        "classification_based_on__classification__variant",
                        "classification_based_on__published_evidence__genome_build__value",
                        "classification_based_on__classification__chgvs_grch37",
                        "classification_based_on__classification__chgvs_grch38",
                    ],
                    renderer=self.render_c_hgvs, client_renderer='TableFormat.hgvs',
                    search=[
                        # "clinvar_allele__allele__clingen_allele__id",  # need a string
                        "classification_based_on__classification__chgvs_grch37",
                        "classification_based_on__classification__chgvs_grch38"
                    ]
            ),
            RichColumn("condition",
                       label="Condition Umbrella",
                       client_renderer='VCTable.condition',
                       search=["condition__display_text"],
                       sort_keys=["condition__sort_text"]
            ),
            RichColumn(name="batches", label="In Batche IDs", renderer=self.batches, client_renderer='renderBatches', orderable=False, search=False, extra_columns=["id"]),
            RichColumn("status", label="Sync Status", client_renderer='renderStatus', sort_keys=["status_sort"], orderable=True, search=False),
            RichColumn("scv", label="SCV", orderable=True),
        ]

    def get_initial_query_params(self, clinvar_key: Optional[str] = None, status: Optional[str] = None) -> QuerySet[ClinVarExport]:
        clinvar_keys = ClinVarKey.clinvar_keys_for_user(self.user)
        if clinvar_key:
            clinvar_keys = clinvar_keys.filter(pk=clinvar_key)

        cve = ClinVarExport.objects.filter(clinvar_allele__clinvar_key__in=clinvar_keys)
        if status:
            statuses = status.split(',')
            cve = cve.filter(status__in=statuses)

        return cve

    def get_initial_queryset(self) -> QuerySet[ClinVarExport]:
        initial_qs = self.get_initial_query_params(
            clinvar_key=self.get_query_param('clinvar_key'),
            status=self.get_query_param('status')
        )
        # add sorting for status so important changes show first
        whens = [
            When(status=ClinVarExportStatus.NEW_SUBMISSION, then=Value(1)),
            When(status=ClinVarExportStatus.CHANGES_PENDING, then=Value(2)),
            When(status=ClinVarExportStatus.UP_TO_DATE, then=Value(3)),
            When(status=ClinVarExportStatus.IN_ERROR, then=Value(4)),
            When(status=ClinVarExportStatus.EXCLUDE, then=Value(5))
        ]
        case = Case(*whens, default=Value(0), output_field=IntegerField())
        initial_qs = initial_qs.annotate(status_sort=case)

        return initial_qs


def clinvar_export_review(request: HttpRequest, clinvar_export_id: int) -> HttpResponseBase:
    clinvar_export: ClinVarExport = ClinVarExport.objects.get(pk=clinvar_export_id)  # fixme get or 404
    clinvar_export.clinvar_allele.clinvar_key.check_user_can_access(request.user)

    if request.method == "POST":
        clinvar_export.scv = request.POST.get("scv") or ""
        clinvar_export.save()
        add_save_message(request, valid=True, name="ClinVarExport")
        return redirect(clinvar_export)

    return render(request, 'classification/clinvar_export.html', context={
        "clinvar_export": clinvar_export,
        "cm": clinvar_export.classification_based_on
    })


def clinvar_export_history(request: HttpRequest, clinvar_export_id: int) -> HttpResponseBase:
    clinvar_export: ClinVarExport = ClinVarExport.objects.get(pk=clinvar_export_id)
    clinvar_export.clinvar_allele.clinvar_key.check_user_can_access(request.user)

    history = clinvar_export.clinvarexportsubmission_set.order_by('-created')

    return render(request, 'classification/clinvar_export_history.html', context={
        'clinvar_export': clinvar_export,
        'history': history
    })


class ClinVarExportSummary(ExportRow):

    def __init__(self, clinvar_export: ClinVarExport):
        self.clinvar_export = clinvar_export

    @property
    def classification(self):
        return self.clinvar_export.classification_based_on

    @lazy
    def genome_build_row(self) -> GenomeBuild:
        if classification := self.classification:
            try:
                return classification.get_genome_build()
            except KeyError:
                pass

    @export_column("ID")
    def id(self):
        return self.clinvar_export.id

    @export_column("Allele URL")
    def allele_id(self) -> str:
        if allele_id := self.classification.classification.allele_id:
            return get_url_from_view_path(reverse('view_allele', kwargs={"allele_id": allele_id}))

    @export_column("URL")
    def url(self):
        return get_url_from_view_path(self.clinvar_export.get_absolute_url())

    @export_column("Genome Build")
    def genome_build(self):
        return str(self.genome_build_row) if self.genome_build_row else None

    @export_column("ClinGenAllele")
    def clingen_allele(self):
        if classification := self.classification:
            if allele := classification.classification.allele:
                return str(allele.clingen_allele)

    @export_column("c.HGVS")
    def c_hgvs(self):
        if genome_build := self.genome_build_row:
            return self.classification.classification.get_c_hgvs(genome_build)

    @export_column("Condition Umbrella")
    def condition_umbrella(self):
        return self.clinvar_export.condition_resolved.as_plain_text

    @export_column("Condition")
    def condition(self):
        if classification := self.classification:
            return classification.condition_text

    @export_column("Interpretation Summary")
    def interpretation_summary(self):
        if classification := self.classification:
            return html_to_text(classification.get(SpecialEKeys.INTERPRETATION_SUMMARY))

    @export_column("Clinical Significance")
    def clinical_significance(self):
        if classification := self.classification:
            return EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))

    @export_column("Sync Status")
    def sync_status(self):
        return self.clinvar_export.get_status_display()

    @export_column("SCV")
    def scv(self):
        return self.clinvar_export.scv

    @lazy
    def latest_submission(self) -> ClinVarExportSubmission:
        return self.clinvar_export.clinvarexportsubmission_set.order_by('-pk').first()

    @export_column("Latest Batch ID")
    def batch_id(self):
        if submission := self.latest_submission:
            return submission.submission_batch_id

    @export_column("Latest Batch Status")
    def batch_status(self):
        if submission := self.latest_submission:
            submission.submission_batch.get_status_display()


def clinvar_export_download(request: HttpRequest, clinvar_key_id: str) -> HttpResponseBase:
    clinvar_key: ClinVarKey = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
    clinvar_key.check_user_can_access(request.user)

    def rows() -> Iterable[str]:
        return ClinVarExportSummary.csv_generator(
            ClinVarExport.objects.filter(clinvar_allele__clinvar_key=clinvar_key).order_by('-id')
        )

    date_str = local_date_string()
    response = StreamingHttpResponse(rows(), content_type='text/csv')
    response['Content-Disposition'] = f'attachment; filename="export_batches_{date_str}.csv"'
    return response


def clinvar_export_batch_download(request: HttpRequest, clinvar_export_batch_id: int) -> HttpResponseBase:
    clinvar_export_batch: ClinVarExportBatch = ClinVarExportBatch.objects.get(pk=clinvar_export_batch_id)
    clinvar_export_batch.clinvar_key.check_user_can_access(request.user)

    # code duplicated from admin, but don't feel this should go into models
    batch_json = clinvar_export_batch.to_json()
    batch_json_str = json.dumps(batch_json)
    response = HttpResponse(batch_json_str, content_type='application/json')
    response['Content-Disposition'] = f'attachment; filename=clinvar_export_preview_{clinvar_export_batch.pk}.json'
    return response


def clinvar_export_summary(request: HttpRequest, clinvar_key_id: Optional[str] = None) -> HttpResponseBase:
    clinvar_key: ClinVarKey
    if not clinvar_key_id:
        if clinvar_key := ClinVarKey.clinvar_keys_for_user(request.user).first():
            return redirect(reverse('clinvar_key_summary', kwargs={'clinvar_key_id': clinvar_key.pk}))
        else:
            # page has special support if no clinvar key is available to the user
            return render(request, 'classification/clinvar_key_summary_none.html')

    clinvar_key = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
    clinvar_key.check_user_can_access(request.user)

    labs = Lab.objects.filter(clinvar_key=clinvar_key).order_by('name')
    dashboard = ClassificationDashboard(user=request.user, labs=labs)

    export_columns = ClinVarExportColumns(request)
    export_batch_columns = ClinVarExportBatchColumns(request)

    return render(request, 'classification/clinvar_key_summary.html', {
        'all_keys': ClinVarKey.clinvar_keys_for_user(request.user),
        'clinvar_key': clinvar_key,
        'labs': labs,
        'missing_condition_count': dashboard.classifications_wout_standard_text,
        'count_records': export_columns.get_initial_query_params(clinvar_key=clinvar_key_id).count(),
        'count_batch': export_batch_columns.get_initial_query_params(clinvar_key=clinvar_key_id).count()
    })


def clinvar_export_detail(request: HttpRequest, clinvar_export_id: int) -> HttpResponseBase:
    clinvar_export = get_object_or_404(ClinVarExport, pk=clinvar_export_id)
    clinvar_export.clinvar_allele.clinvar_key.check_user_can_access(request.user)

    interpretation_summary: Optional[str] = None
    if cm := clinvar_export.classification_based_on:
        if interpret_summary_html := cm.get(SpecialEKeys.INTERPRETATION_SUMMARY):
            interpretation_summary = html_to_text(interpret_summary_html)

    return render(request, 'classification/clinvar_export_detail.html', {
        'clinvar_export': clinvar_export,
        'classification': clinvar_export.classification_based_on,
        'interpretation_summary': interpretation_summary
    })


class ClinVarMatchView(View):

    def get(self, request, **kwargs):
        clinvar_key_id = kwargs.get('clinvar_key_id')
        clinvar_key: ClinVarKey
        if not clinvar_key_id:
            if clinvar_key := ClinVarKey.clinvar_keys_for_user(request.user).first():
                return redirect(reverse('clinvar_match', kwargs={'clinvar_key_id': clinvar_key.pk}))
            else:
                # page has special support if no clinvar key is available to the user
                return render(request, 'classification/clinvar_key_summary_none.html')

        clinvar_key = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
        clinvar_key.check_user_can_access(request.user)

        return render(request, 'classification/clinvar_match.html', {
            'all_keys': ClinVarKey.clinvar_keys_for_user(request.user),
            'clinvar_key': clinvar_key
        })

    def post(self, request, **kwargs):
        clinvar_key_id = kwargs.get('clinvar_key_id')
        clinvar_key: ClinVarKey = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
        clinvar_key.check_user_can_access(request.user)

        file_obj = io_string = io.StringIO(request.FILES.get('file').read().decode("utf-8"))
        clinvar_legacy_rows = ClinVarLegacyRow.load_file(file_obj, clinvar_key)

        return render(request, 'classification/clinvar_match.html', {
            'all_keys': ClinVarKey.clinvar_keys_for_user(request.user),
            'clinvar_key': clinvar_key,
            'rows': clinvar_legacy_rows
        })


def clinvar_match_detail(request, clinvar_key_id: str):
    clinvar_key: ClinVarKey = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
    clinvar_key.check_user_can_access(request.user)

    # TODO should this be a post?
    data_str = request.GET.get('data_str')

    return render(request, 'classification/clinvar_match_detail.html', {
        'matches': ClinVarLegacyRow.from_data_str(clinvar_key, data_str).find_variant_grid_allele()
    })