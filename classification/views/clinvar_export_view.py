from typing import Dict, Any

from django.shortcuts import render
from guardian.shortcuts import get_objects_for_user
from gunicorn.config import User
from htmlmin.decorators import not_minified_response

from classification.models import ClinVarExport, EvidenceKeyMap, ConditionTextMatch, ClassificationModification
from snpdb.views.datatable_view import DatatableConfig, RichColumn


class ClinVarExportColumns(DatatableConfig):

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

    def __init__(self, request):
        super().__init__(request)

        evidence_keys = EvidenceKeyMap.cached()

        self.rich_columns = [
            RichColumn(key="classification_based_on__created", label='ClinVar Variant', orderable=True,
                       sort_keys=["classification_based_on__classification__chgvs_grch38"], extra_columns=[
                "id",
                "classification_based_on__created",
                "classification_based_on__classification_id",
                "classification_based_on__published_evidence__genome_build__value",
                "classification_based_on__classification__chgvs_grch37",
                "classification_based_on__classification__chgvs_grch38",
            ], renderer=self.render_classification_link, client_renderer='renderId'),

            RichColumn('lab__name', name='Lab', orderable=True),
            RichColumn('gene_symbol', name='Gene Symbol', orderable=True),

            # this busy ups the table a little too much
            # evidence_keys.get(SpecialEKeys.MODE_OF_INHERITANCE).rich_column(prefix="classification_based_on__published_evidence"),

            # RichColumn('condition_text_normal', label='Condition Text', orderable=True),
            # RichColumn('condition_xrefs', label='Terms', orderable=True, client_renderer='ontologyList', renderer=self.render_aliases),
            RichColumn('submit_when_possible', name='Auto-Submit Enabled', orderable=True, client_renderer='TableFormat.boolean.bind(null, "standard")')
        ]

    def get_initial_queryset(self):
        return get_objects_for_user(self.user, ClinVarExport.get_read_perm(), klass=ClinVarExport, accept_global_perms=True)


def clinvar_exports_view(request):
    return render(request, 'classification/clinvar_exports.html', context={
        'datatable_config': ClinVarExportColumns(request)
    })


@not_minified_response
def clinvar_export_review_view(request, pk):
    user: User = request.user
    clinvar_export: ClinVarExport = ClinVarExport.objects.get(pk=pk)
    condition_text_match = ConditionTextMatch.objects.filter(classification=clinvar_export.classification_based_on.classification).first()
    cm: ClassificationModification = clinvar_export.classification_based_on

    return render(request, 'classification/clinvar_export.html', context={
        'clinvar_export': clinvar_export,
        'condition_text_match': condition_text_match,
        'cm': cm
    })
