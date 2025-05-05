import io
from functools import cached_property

from django.contrib import messages
from django.db.models import QuerySet
from django.shortcuts import redirect, get_object_or_404, render
from django.urls import reverse
from django.views import View

from classification.models.clinvar_legacy import ClinVarLegacy
from classification.services.clinvar_legacy_services import ClinVarLegacyService
from library.django_utils import RequireSuperUserView
from snpdb.models import ClinVarKey
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder, DC


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


class ClinVarLegacyColumns(DatatableConfig[ClinVarLegacy]):

    def __init__(self, request: any):
        super().__init__(request)
        self.rich_columns = [
            RichColumn("scv", default_sort=SortOrder.ASC),
            RichColumn("clinvar_allele_id", label="ClinVar Allele"),
        ]

    @cached_property
    def get_clinvar_key(self):
        return ClinVarKey.objects.filter(id=self.get_query_param('clinvar_key_id')).first()

    def get_initial_queryset(self) -> QuerySet[DC]:
        return ClinVarLegacy.objects.filter(clinvar_key=self.get_clinvar_key)

