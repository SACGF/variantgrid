from typing import Dict, Any

from django.core.exceptions import PermissionDenied
from django.shortcuts import render
from gunicorn.config import User
from htmlmin.decorators import not_minified_response

from classification.models import EvidenceKeyMap, ClinVarExportRecord, ClinVarAllele
from library.cache import timed_cache
from snpdb.models import Allele
from snpdb.views.datatable_view import DatatableConfig, RichColumn

@timed_cache(size_limit=30, ttl=60)
def allele_for(allele_id: int):
    return Allele.objects.select_related('clingen_allele').get(pk=allele_id)

class ClinVarExportRecordColumns(DatatableConfig):

    def render_classification_link(self, row: Dict[str, Any]):
        created = row["classification_based_on__created"]
        c_id = row["classification_based_on__classification_id"]

        link_id = f"{c_id}.{created.timestamp()}"

        genome_build = row["classification_based_on__published_evidence__genome_build__value"]
        c_hgvs = None
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

    def __init__(self, request):
        super().__init__(request)

        evidence_keys = EvidenceKeyMap.cached()

        self.rich_columns = [
            # FIXME this is all based on the previous stuff
            RichColumn("clinvar_allele__clinvar_key", name="ClinVar Key"),
            RichColumn("clinvar_allele__allele_id", renderer=self.render_allele, name="Allele"),
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
        return ClinVarExportRecord.objects.filter(clinvar_allele__clinvar_key__in=ClinVarAllele.clinvar_keys_for_user(self.user))
        # return get_objects_for_user(self.user, ClinVarExport.get_read_perm(), klass=ClinVarExport, accept_global_perms=True)


def clinvar_exports_view(request):
    return render(request, 'classification/clinvar_exports.html', context={
        'datatable_config': ClinVarExportRecordColumns(request)
    })


@not_minified_response
def clinvar_export_review_view(request, pk):
    clinvar_export: ClinVarExportRecord = ClinVarExportRecord.objects.get(pk=pk)  # fixme get or 404

    user: User = request.user
    if not user.is_superuser:
        allowed_clinvar_keys = ClinVarAllele.clinvar_keys_for_user(user)
        if not allowed_clinvar_keys.filter(pk=clinvar_export.clinvar_allele.clinvar_key).exists():
            raise PermissionDenied("User does not belong to a lab that uses the submission key")

    return render(request, 'classification/clinvar_export.html', context={
        "clinvar_export": clinvar_export,
        "cm": clinvar_export.classification_based_on
    })