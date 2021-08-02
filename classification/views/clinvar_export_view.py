from typing import Dict, Any, Optional

from django.contrib import messages
from django.http import HttpResponse
from django.shortcuts import render, redirect, get_object_or_404
from django.urls import reverse
from htmlmin.decorators import not_minified_response

from classification.enums import ShareLevel, SpecialEKeys
from classification.models import ClinVarExport, ClinVarExportSubmissionBatch, ClinVarExportStatus, ClinVarAllele, \
    ClassificationModification, ClinVarExportSubmissionbatchStatus, Classification
from library.cache import timed_cache
from library.django_utils import add_save_message
from library.utils import html_to_text
from snpdb.models import Allele, ClinVarKey, Lab
from snpdb.views.datatable_view import DatatableConfig, RichColumn
import json


@timed_cache(size_limit=30, ttl=60)
def allele_for(allele_id: int) -> Allele:
    return Allele.objects.select_related('clingen_allele').get(pk=allele_id)


class ClinVarExportAlleleColumns(DatatableConfig):

    """
        allele = models.ForeignKey(Allele, on_delete=models.CASCADE)
    clinvar_key = models.ForeignKey(ClinVarKey, null=True, blank=True, on_delete=models.CASCADE)

    classificatiocns_missing_condition = models.IntegerField(default=0)
    submissions_valid = models.IntegerField(default=0)
    submissions_invalid = models.IntegerField(default=0)
    last_evaluated = models.DateTimeField(default=now)
    """

    def render_allele(self, row: Dict[str, Any]) -> str:
        allele_id = row.get('allele')
        allele = allele_for(allele_id)
        return f"{allele:CA}"

    def __init__(self, request):
        super().__init__(request)

        self.expand_client_renderer = "TableFormat.expandAjax.bind(null, 'clinvar_export_alleles_datatable_expand', 'id')";
        self.rich_columns = [
            RichColumn("id", label="ID", orderable=True),
            RichColumn("clinvar_key", label="ClinVar Key", orderable=True),
            RichColumn("allele", renderer=self.render_allele, orderable=True),
            RichColumn("last_evaluated", client_renderer='TableFormat.timeAgo', orderable=True),
            RichColumn("submissions_valid", client_renderer='TableFormat.severeNumber.bind(null, "success")', label="Submissions Valid", orderable=True),
            RichColumn("submissions_invalid", client_renderer='TableFormat.severeNumber.bind(null, "danger")', label="Submission w Errors", orderable=True),
            RichColumn("classifications_missing_condition", client_renderer='TableFormat.severeNumber.bind(null, "danger")', label="Classifications w/out Standard Condition", orderable=True),
        ]

    def get_initial_queryset(self):
        return ClinVarAllele.objects.filter(clinvar_key__in=ClinVarKey.clinvar_keys_for_user(self.user))


class ClinVarExportBatchColumns(DatatableConfig):

    def render_status(self, row: Dict[str, Any]):
        return ClinVarExportSubmissionbatchStatus(row['status']).label

    def __init__(self, request):
        super().__init__(request)

        """
        clinvar_key = models.ForeignKey(ClinVarKey, on_delete=models.PROTECT)
        submission_version = models.IntegerField()
        status = models.CharField(max_length=1, choices=ClinVarExportSubmissionbatchStatus.choices, default=ClinVarExportSubmissionbatchStatus.AWAITING_UPLOAD)
        submitted_json = models.JSONField(null=True, blank=True)  # leave this as blank until we've actually uploaded data
        response_json = models.JSONField(null=True, blank=True)  # what did ClinVar return
        """

        self.expand_client_renderer = "TableFormat.expandAjax.bind(null, 'clinvar_export_batches_datatable_expand', 'id')";
        self.rich_columns = [
            # TODO, toggle enabled on screens where it makes sense
            RichColumn("id", orderable=True),
            RichColumn("clinvar_key", name="ClinVar Key", orderable=True, enabled=False),
            RichColumn("created", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn("status", renderer=self.render_status, orderable=True)
        ]

    def get_initial_query_params(self, clinvar_key: Optional[str] = None):
        clinvar_keys = ClinVarKey.clinvar_keys_for_user(self.user)
        if clinvar_key:
            clinvar_keys = clinvar_keys.filter(pk=clinvar_key)
        cve = ClinVarExportSubmissionBatch.objects.filter(clinvar_key__in=clinvar_keys)
        return cve

    def get_initial_queryset(self):
        return self.get_initial_query_params(
            clinvar_key=self.get_query_param('clinvar_key')
        )


def clinvar_export_batches_datatable_expand(request, pk: int):
    clinvar_export_batch: ClinVarExportSubmissionBatch = ClinVarExportSubmissionBatch.objects.get(pk=pk)
    clinvar_export_batch.clinvar_key.check_user_can_access(request.user)
    submissions = clinvar_export_batch.clinvarexportsubmission_set.order_by('created')

    return render(request, 'classification/clinvar_export_batch_expand.html', {
        "batch": clinvar_export_batch,
        "submissions": submissions
    })


class ClinVarExportRecordColumns(DatatableConfig):

    def render_classification_link(self, row: Dict[str, Any]):
        created = row["classification_based_on__created"]
        c_id = row["classification_based_on__classification_id"]

        link_id = f"{c_id}.{created.timestamp()}"

        genome_build = row["classification_based_on__published_evidence__genome_build__value"]
        c_hgvs: str
        if "h37" in genome_build:
            c_hgvs = row["classification_based_on__classification__chgvs_grch37"]
        else:
            c_hgvs = row["classification_based_on__classification__chgvs_grch38"]

        return {
            "id": row["id"],
            "genome_build": genome_build,
            "c_hgvs": c_hgvs,
            "cm_id": link_id
        }

    def render_allele(self, row: Dict[str, Any]) -> str:
        allele_id = row.get('clinvar_allele__allele_id')
        allele = allele_for(allele_id)
        return f"{allele:CA}"

    def render_status(self, row: Dict[str, Any]):
        return ClinVarExportStatus(row['status']).label

    def __init__(self, request):
        super().__init__(request)

        self.expand_client_renderer = "TableFormat.expandAjax.bind(null, 'clinvar_export_detail', 'id')";
        self.rich_columns = [
            # TODO, toggle enabled on screens where it makes sense
            RichColumn("id", orderable=True),
            RichColumn("clinvar_allele__clinvar_key", name="ClinVar Key", orderable=True, enabled=False),
            RichColumn("clinvar_allele__allele_id", renderer=self.render_allele, name="allele", label="Allele", orderable=True),
            RichColumn("status", renderer=self.render_status, orderable=True),
            RichColumn("scv", label="SCV", orderable=True),
            RichColumn(key="classification_based_on__created", name="classification", label='ClinVar Variant',
                       sort_keys=["classification_based_on__classification__chgvs_grch38"], extra_columns=[
                        "id",
                        "classification_based_on__created",
                        "classification_based_on__classification_id",
                        "classification_based_on__published_evidence__genome_build__value",
                        "classification_based_on__classification__chgvs_grch37",
                        "classification_based_on__classification__chgvs_grch38",
                        ], renderer=self.render_classification_link, client_renderer='renderId'),
            RichColumn("condition", name="condition", client_renderer='VCTable.condition'),

            # RichColumn('lab__name', name='Lab', orderable=True),

            # this busy ups the table a little too much
            # evidence_keys.get(SpecialEKeys.MODE_OF_INHERITANCE).rich_column(prefix="classification_based_on__published_evidence"),

            # RichColumn('condition_text_normal', label='Condition Text', orderable=True),
            # RichColumn('condition_xrefs', label='Terms', orderable=True, client_renderer='ontologyList', renderer=self.render_aliases),
            # RichColumn('submit_when_possible', name='Auto-Submit Enabled', orderable=True, client_renderer='TableFormat.boolean.bind(null, "standard")')
        ]

    def get_initial_query_params(self, clinvar_key: Optional[str] = None, status: Optional[str] = None):
        clinvar_keys = ClinVarKey.clinvar_keys_for_user(self.user)
        if clinvar_key:
            clinvar_keys = clinvar_keys.filter(pk=clinvar_key)

        cve = ClinVarExport.objects.filter(clinvar_allele__clinvar_key__in=clinvar_keys)
        if status:
            statuses = status.split(',')
            cve = cve.filter(status__in=statuses)

        return cve

    def get_initial_queryset(self):
        return self.get_initial_query_params(
            clinvar_key=self.get_query_param('clinvar_key'),
            status=self.get_query_param('status')
        )


def clinvar_export_alleles_view(request):
    return render(request, 'classification/clinvar_export_alleles.html', context={
        'datatable_config': ClinVarExportAlleleColumns(request)
    })


def clinvar_export_allele_datatable_expand_view(request, pk: int):
    clinvar_allele: ClinVarAllele = ClinVarAllele.objects.get(pk=pk)
    clinvar_allele.clinvar_key.check_user_can_access(request.user)

    submissions = ClinVarExport.objects.filter(clinvar_allele=clinvar_allele).order_by('status')
    labs = clinvar_allele.clinvar_key.lab_set.all()
    missing_conditions = ClassificationModification.latest_for_user(user=request.user, allele=clinvar_allele.allele, published=True, shared_only=True).filter(classification__condition_resolution__isnull=True, classification__lab__in=labs)
    reduced_missing_conditions = clinvar_allele.classifications_missing_condition - len(missing_conditions)

    return render(request, 'classification/clinvar_export_allele_expand.html', {
        "allele": clinvar_allele.allele,
        "submissions": submissions,
        "missing_conditions": missing_conditions,
        "reduced_missing_conditions": reduced_missing_conditions
    })


def clinvar_exports_view(request):
    return render(request, 'classification/clinvar_exports.html', context={
        'datatable_config': ClinVarExportRecordColumns(request)
    })


@not_minified_response
def clinvar_export_review_view(request, pk):
    clinvar_export: ClinVarExport = ClinVarExport.objects.get(pk=pk)  # fixme get or 404
    clinvar_export.clinvar_allele.clinvar_key.check_user_can_access(request.user)

    if request.method == "POST":
        clinvar_export.scv = request.POST.get("scv") or None
        clinvar_export.save()
        add_save_message(request, valid=True, name="ClinVarExport")
        return redirect(clinvar_export)

    return render(request, 'classification/clinvar_export.html', context={
        "clinvar_export": clinvar_export,
        "cm": clinvar_export.classification_based_on
    })


@not_minified_response
def clinvar_export_history_view(request, pk):
    clinvar_export: ClinVarExport = ClinVarExport.objects.get(pk=pk)
    clinvar_export.clinvar_allele.clinvar_key.check_user_can_access(request.user)

    history = clinvar_export.clinvarexportsubmission_set.order_by('-created')

    return render(request, 'classification/clinvar_export_history.html', context={
        'clinvar_export': clinvar_export,
        'history': history
    })


def clinvar_export_batch_view(request, pk):
    clinvar_export_batch: ClinVarExportSubmissionBatch = ClinVarExportSubmissionBatch.objects.get(pk=pk)
    clinvar_export_batch.clinvar_key.check_user_can_access(request.user)

    return render(request, 'classification/clinvar_export_batch.html', context={
        'batch': clinvar_export_batch,
        'submissions': clinvar_export_batch.clinvarexportsubmission_set.order_by('created')
    })


def clinvar_export_batch_download(request, pk):
    clinvar_export_batch: ClinVarExportSubmissionBatch = ClinVarExportSubmissionBatch.objects.get(pk=pk)
    clinvar_export_batch.clinvar_key.check_user_can_access(request.user)

    # code duplicated from admin, but don't feel this should go into models
    batch_json = clinvar_export_batch.to_json()
    batch_json_str = json.dumps(batch_json)
    response = HttpResponse(batch_json_str, content_type='application/json')
    response['Content-Disposition'] = f'attachment; filename=clinvar_export_preview_{clinvar_export_batch.pk}.json'
    return response


def clinvar_export_summary(request, pk: Optional[str] = None):
    clinvar_key: ClinVarKey
    if not pk:
        clinvar_key = ClinVarKey.clinvar_keys_for_user(request.user).first()
        return redirect(reverse('clinvar_key_summary', kwargs={'pk':clinvar_key.pk}))

    clinvar_key = get_object_or_404(ClinVarKey, pk=pk)
    clinvar_key.check_user_can_access(request.user)

    labs = Lab.objects.filter(clinvar_key=clinvar_key).order_by('name')
    missing_condition = Classification.objects.filter(lab__in=labs, share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS, condition_resolution__isnull=True)

    export_columns = ClinVarExportRecordColumns(request)
    export_batch_columns = ClinVarExportBatchColumns(request)

    return render(request, 'classification/clinvar_key_summary.html', {
        'all_keys': ClinVarKey.clinvar_keys_for_user(request.user),
        'clinvar_key': clinvar_key,
        'labs': labs,
        'missing_condition_count': missing_condition.count(),
        'export_columns': export_columns,
        'export_batch_columns': export_batch_columns,
        'count_records': export_columns.get_initial_query_params(clinvar_key=pk).count(),
        'count_batch': export_batch_columns.get_initial_query_params(clinvar_key=pk).count()
    })


def clinvar_export_expand(request, pk: Optional[str] = None):
    clinvar_export = get_object_or_404(ClinVarExport, pk=pk)
    clinvar_export.clinvar_allele.clinvar_key.check_user_can_access(request.user)

    interpretation_summary: Optional[str] = None
    if cm := clinvar_export.classification_based_on:
        if interpret_summary_html := clinvar_export.classification_based_on.get(SpecialEKeys.INTERPRETATION_SUMMARY):
            interpretation_summary = html_to_text(interpret_summary_html)

    return render(request, 'classification/clinvar_export_expand.html', {
        'clinvar_export': clinvar_export,
        'classification': clinvar_export.classification_based_on,
        'interpretation_summary': interpretation_summary
    })