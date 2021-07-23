from typing import Dict, Any

from django.http import HttpResponse
from django.shortcuts import render
from htmlmin.decorators import not_minified_response

from classification.models import ClinVarExport, ClinVarExportSubmissionBatch, ClinVarStatus
from library.cache import timed_cache
from snpdb.models import Allele, ClinVarKey
from snpdb.views.datatable_view import DatatableConfig, RichColumn
import json


@timed_cache(size_limit=30, ttl=60)
def allele_for(allele_id: int):
    return Allele.objects.select_related('clingen_allele').get(pk=allele_id)


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

    def render_allele(self, row: Dict[str, Any]):
        allele_id = row.get('clinvar_allele__allele_id')
        allele = allele_for(allele_id)
        return str(allele)

    def render_status(self, row: Dict[str, Any]):
        return ClinVarStatus(row['status']).label

    def __init__(self, request):
        super().__init__(request)

        self.rich_columns = [
            # FIXME this is all based on the previous stuff
            RichColumn("clinvar_allele__clinvar_key", name="ClinVar Key"),
            RichColumn("clinvar_allele__allele_id", renderer=self.render_allele, name="Allele"),
            RichColumn("status", renderer=self.render_status),
            RichColumn(key="classification_based_on__created", label='ClinVar Variant', orderable=True,
                       sort_keys=["classification_based_on__classification__chgvs_grch38"], extra_columns=[
                        "id",
                        "classification_based_on__created",
                        "classification_based_on__classification_id",
                        "classification_based_on__published_evidence__genome_build__value",
                        "classification_based_on__classification__chgvs_grch37",
                        "classification_based_on__classification__chgvs_grch38",
                        ], renderer=self.render_classification_link, client_renderer='renderId'),
            RichColumn("condition", name="Condition", client_renderer='VCTable.condition'),

            # RichColumn('lab__name', name='Lab', orderable=True),

            # this busy ups the table a little too much
            # evidence_keys.get(SpecialEKeys.MODE_OF_INHERITANCE).rich_column(prefix="classification_based_on__published_evidence"),

            # RichColumn('condition_text_normal', label='Condition Text', orderable=True),
            # RichColumn('condition_xrefs', label='Terms', orderable=True, client_renderer='ontologyList', renderer=self.render_aliases),
            # RichColumn('submit_when_possible', name='Auto-Submit Enabled', orderable=True, client_renderer='TableFormat.boolean.bind(null, "standard")')
        ]

    def get_initial_queryset(self):
        return ClinVarExport.objects.filter(clinvar_allele__clinvar_key__in=ClinVarKey.clinvar_keys_for_user(self.user))
        # return get_objects_for_user(self.user, ClinVarExport.get_read_perm(), klass=ClinVarExport, accept_global_perms=True)


def clinvar_exports_view(request):
    return render(request, 'classification/clinvar_exports.html', context={
        'datatable_config': ClinVarExportRecordColumns(request)
    })


@not_minified_response
def clinvar_export_review_view(request, pk):
    clinvar_export: ClinVarExport = ClinVarExport.objects.get(pk=pk)  # fixme get or 404
    clinvar_export.clinvar_allele.clinvar_key.check_user_can_access(request.user)

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
