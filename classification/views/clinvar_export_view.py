import json
from typing import Dict, Any, Optional, Iterable

from django.conf import settings
from django.db.models import QuerySet, When, Value, Case, IntegerField, Count
from django.http import HttpResponse, StreamingHttpResponse, HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render, redirect, get_object_or_404
from django.urls import reverse
from django.utils.timezone import now
from lazy import lazy
from pytz import timezone

from classification.enums import SpecialEKeys
from classification.models import ClinVarExport, ClinVarExportBatch, ClinVarExportBatchStatus, \
    EvidenceKeyMap, ClinVarExportStatus
from classification.views.classification_dashboard_view import ClassificationDashboard
from genes.hgvs import CHGVS
from library.cache import timed_cache
from library.django_utils import add_save_message, get_url_from_view_path
from library.utils import html_to_text, export_column, ExportRow
from snpdb.models import ClinVarKey, Lab, Allele
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder
from uicore.json.json_types import JsonDataType


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
    def genome_build_row(self):
        if classification := self.classification:
            try:
                return classification.get_genome_build()
            except KeyError:
                pass

    @export_column("ID")
    def id(self):
        return self.clinvar_export.id

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


def clinvar_export_download(request: HttpRequest, clinvar_key_id: str) -> HttpResponseBase:
    clinvar_key: ClinVarKey = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
    clinvar_key.check_user_can_access(request.user)

    def rows() -> Iterable[str]:
        return ClinVarExportSummary.csv_generator(
            ClinVarExport.objects.filter(clinvar_allele__clinvar_key=clinvar_key).order_by('-id')
        )

    date_str = now().astimezone(tz=timezone(settings.TIME_ZONE)).strftime("%Y-%m-%d")
    response = StreamingHttpResponse(rows(), content_type='text/csv')
    response['Content-Disposition'] = f'attachment; filename="export_batches_{date_str}.csv"'  # TODO  add in date
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
            return redirect(reverse('clinvar_key_summary', kwargs={'pk': clinvar_key.pk}))
        else:
            # page has special support if no clinvar key is available to the user
            return render(request, 'classification/clinvar_key_summary_none.html')

    clinvar_key = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
    clinvar_key.check_user_can_access(request.user)

    labs = Lab.objects.filter(clinvar_key=clinvar_key).order_by('name')
    dashbaord = ClassificationDashboard(user=request.user, labs=labs)

    export_columns = ClinVarExportColumns(request)
    export_batch_columns = ClinVarExportBatchColumns(request)

    return render(request, 'classification/clinvar_key_summary.html', {
        'all_keys': ClinVarKey.clinvar_keys_for_user(request.user),
        'clinvar_key': clinvar_key,
        'labs': labs,
        'missing_condition_count': dashbaord.classifications_wout_standard_text,
        'count_records': export_columns.get_initial_query_params(clinvar_key=pk).count(),
        'count_batch': export_batch_columns.get_initial_query_params(clinvar_key=pk).count()
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
