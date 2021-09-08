"""
@see https://django-autocomplete-light.readthedocs.io/en/master/
"""

import operator
from functools import reduce

from django.db.models.query_utils import Q
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie

from library.constants import MINUTE_SECS
from library.django_utils.autocomplete_utils import AutocompleteView
from patients.models import Patient, Clinician, Specimen, ExternalPK


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class PatientAutocompleteView(AutocompleteView):
    fields = ['last_name', 'first_name']

    def get_user_queryset(self, user):
        return Patient.filter_for_user(user)


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class SpecimenAutocompleteView(AutocompleteView):
    fields = ['reference_id']

    def get_user_queryset(self, user):
        patients_qs = Patient.filter_for_user(user)
        patient = self.forwarded.get('patient', None)
        if patient:
            # print(f"Filtering for forwarded patient: {patient})")
            patients_qs = patients_qs.filter(pk=patient)

        return Specimen.objects.filter(patient__in=patients_qs)


@method_decorator(cache_page(MINUTE_SECS), name='dispatch')
class ClinicianAutocompleteView(AutocompleteView):
    fields = ['last_name', 'first_name']

    def get_user_queryset(self, user):
        return Clinician.objects.all()


@method_decorator(cache_page(30), name='dispatch')
class ExternalPKAutocompleteView(AutocompleteView):
    fields = ['code']

    def get_user_queryset(self, user):
        external_type = self.forwarded.get('external_type', None)

        q_list = []
        if external_type:
            q_list.append(Q(external_type=external_type))

        qs = ExternalPK.objects.all()
        if q_list:
            q = reduce(operator.and_, q_list)
            qs = qs.filter(q)
        return qs
