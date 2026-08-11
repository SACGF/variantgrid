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
from patients.models import Clinician, Extraction, ExternalPK, Patient, Specimen


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class PatientAutocompleteView(AutocompleteView):
    fields = ['last_name', 'first_name']

    def get_user_queryset(self, user):
        return Patient.filter_for_user(user)


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class SpecimenAutocompleteView(AutocompleteView):
    """ Narrows on whichever of Patient -> Specimen -> Extraction a form has already set """
    fields = ['reference_id']

    def get_user_queryset(self, user):
        qs = Specimen.objects.filter(patient__in=Patient.filter_for_user(user))
        if patient := self.forwarded.get('patient'):
            qs = qs.filter(patient=patient)
        if extraction := self.forwarded.get('extraction'):
            qs = qs.filter(extraction=extraction)
        return qs


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class ExtractionAutocompleteView(AutocompleteView):
    """ Narrows on whichever of Patient -> Specimen a form has already set """
    # An extraction may be referred to by its own reference (eg a container suffix) or its specimen's
    fields = ['reference_id', 'specimen__reference_id']

    def get_user_queryset(self, user):
        qs = Extraction.objects.filter(specimen__patient__in=Patient.filter_for_user(user))
        if patient := self.forwarded.get('patient'):
            qs = qs.filter(specimen__patient=patient)
        if specimen := self.forwarded.get('specimen'):
            qs = qs.filter(specimen=specimen)
        return qs


@method_decorator(cache_page(MINUTE_SECS), name='dispatch')
class ClinicianAutocompleteView(AutocompleteView):
    fields = ['last_name', 'first_name']

    def get_user_queryset(self, user):
        return Clinician.objects.all()


@method_decorator(cache_page(30), name='dispatch')
class ExternalPKAutocompleteView(AutocompleteView):
    fields = ['code']

    def get_user_queryset(self, _user):
        """ This is just a PK with no identifiable info, so doesn't require
            any per-user filtering """
        external_type = self.forwarded.get('external_type', None)
        q_list = []
        if external_type:
            q_list.append(Q(external_type=external_type))

        qs = ExternalPK.objects.all()
        if q_list:
            q = reduce(operator.and_, q_list)
            qs = qs.filter(q)
        return qs
