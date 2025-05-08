import io
import re
from functools import cached_property
from django.contrib import messages
from django.db.models import QuerySet
from django.shortcuts import redirect, get_object_or_404, render
from django.urls import reverse
from django.views import View
from classification.models import ConditionResolved, ClinVarExport
from classification.models.clinvar_legacy import ClinVarLegacy
from classification.services.clinvar_legacy_services import ClinVarLegacyService
from genes.hgvs import CHGVS
from library.django_utils import RequireSuperUserView
from library.utils import JsonDataType
from library.utils.django_utils import render_ajax_view
from ontology.models import OntologyTerm
from snpdb.models import ClinVarKey
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder, DC, CellData, DatatableConfigQuerySetMode


class ClinVarLegacyView(RequireSuperUserView, View):

    def get(self, request, **kwargs):
        clinvar_key_id = kwargs.get('clinvar_key_id')
        clinvar_key: ClinVarKey
        if not clinvar_key_id:
            if clinvar_key := ClinVarKey.clinvar_keys_for_user(request.user).first():
                return redirect(reverse('clinvar_legacy', kwargs={'clinvar_key_id': clinvar_key.pk}))
            else:
                # page has special support if no clinvar key is available to the user
                return render(request, 'classification/clinvar_key_summary_none.html')

        clinvar_key = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
        clinvar_key.check_user_can_access(request.user)

        return render(request, 'classification/clinvar_legacy.html', {
            'all_keys': ClinVarKey.clinvar_keys_for_user(request.user),
            'clinvar_key': clinvar_key
        })

    def post(self, request, **kwargs):
        clinvar_key_id = kwargs.get('clinvar_key_id')
        clinvar_key: ClinVarKey = get_object_or_404(ClinVarKey, pk=clinvar_key_id)
        clinvar_key.check_user_can_access(request.user)

        try:
            file_obj = io.StringIO(request.FILES.get('file').read().decode("utf-8"))
            filename = request.FILES.get('file').name
            delimiter = "\t" if filename.endswith(".tsv") else ","
            ClinVarLegacyService(clinvar_key=clinvar_key).load_file(file_obj, delimiter=delimiter)
        except AttributeError as ve:
            messages.error(request, f"File is missing or in the wrong format ({ve})")

        return render(request, 'classification/clinvar_legacy.html', {
            'all_keys': ClinVarKey.clinvar_keys_for_user(request.user),
            'clinvar_key': clinvar_key
        })


def view_clinvar_legacy_detail(request, scv: str):
    clinvar_legacy = get_object_or_404(ClinVarLegacy, scv=scv)
    clinvar_legacy.clinvar_key.check_user_can_access(request.user)
    matches = ClinVarLegacyService.find_matches(clinvar_legacy)
    return render_ajax_view(request, 'classification/clinvar_legacy_detail.html', {
        "clinvar_legacy": clinvar_legacy,
        "matches": matches,
    })


class ClinVarLegacyColumns(DatatableConfig[ClinVarLegacy]):

    def render_hgvs(self, row: CellData[ClinVarLegacy]) -> JsonDataType:
        return row.obj.hgvs_obj.to_json()

    def render_condition(self, row: CellData[ClinVarLegacy]) -> JsonDataType:
        if value := row.obj.your_condition_identifier:
            stub = OntologyTerm.get_or_stub(value)
            return ConditionResolved(terms=[stub]).to_json()
        return None

    def render_classification(self, row: CellData[ClinVarLegacy]) -> JsonDataType:
        return row.obj.classification_code

    def render_last_evaluated(self, row: CellData[ClinVarLegacy]) -> JsonDataType:
        return row.obj.date_last_evaluated

    def render_scv(self, row: CellData[ClinVarLegacy]) -> JsonDataType:
        return row.obj.scv
        # scv = row.obj.scv
        # if existing := ClinVarExport.objects.filter(scv=scv).first():
        #     return {"url": existing.get_absolute_url(), "text": existing.scv}
        # return {"text": scv}

    def render_allele(self, row: CellData[ClinVarLegacy]) -> JsonDataType:
        if allele := row.obj.allele:
            return {"url": allele.get_absolute_url(), "text": f"A{allele.pk}"}
        else:
            return None


    def __init__(self, request: any):
        super().__init__(request)
        self.expand_client_renderer = DatatableConfig._row_expand_ajax('clinvar_legacy_detail', id_field="scv", expected_height=120)
        self.server_calculate_mode = DatatableConfigQuerySetMode.OBJECTS
        self.rich_columns = [
            RichColumn(name="scv", label="SCV", renderer=self.render_scv, default_sort=SortOrder.ASC, sort_keys=["scv"]),
            RichColumn("allele", label="Allele", renderer=self.render_allele, client_renderer="TableFormat.linkUrl"),
            RichColumn("preferred_variant_name", label="HGVS", renderer=self.render_hgvs, client_renderer="VCTable.hgvs"),
            RichColumn("your_condition_identifier", label="Condition", renderer=self.render_condition, client_renderer="VCTable.condition"),
            RichColumn("clinical_significance", label="Classification", renderer=self.render_classification, client_renderer="VCTable.classification"),
            RichColumn("date_last_evaluated", label="Date Last Evaluated", client_renderer="TableFormat.timestamp"),
        ]

    @cached_property
    def get_clinvar_key(self):
        return ClinVarKey.objects.filter(id=self.get_query_param('clinvar_key_id')).first()

    def get_initial_queryset(self) -> QuerySet[DC]:
        return ClinVarLegacy.objects.filter(clinvar_key=self.get_clinvar_key)

