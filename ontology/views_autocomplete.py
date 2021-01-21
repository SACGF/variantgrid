import abc

from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page

from library.constants import HOUR_SECS
from library.django_utils.autocomplete_utils import AutocompleteView
from ontology.models import OntologyTerm, OntologyService


class AbstractOntologyTermAutocompleteView(abc.ABC, AutocompleteView):
    fields = ['name']

    @abc.abstractmethod
    def _get_ontology_service(self):
        pass

    def get_user_queryset(self, user):
        qs = OntologyTerm.objects.all()
        if ontology_service := self._get_ontology_service():
            qs = qs.filter(ontology_service=ontology_service)
        return qs


@method_decorator(cache_page(HOUR_SECS), name='dispatch')
class OntologyTermAutocompleteView(AbstractOntologyTermAutocompleteView):
    def _get_ontology_service(self):
        # Passed ontology_service in forward
        return self.forwarded.get('ontology_service')


@method_decorator(cache_page(HOUR_SECS), name='dispatch')
class HPOAutocompleteView(AbstractOntologyTermAutocompleteView):
    def _get_ontology_service(self):
        return OntologyService.HPO


@method_decorator(cache_page(HOUR_SECS), name='dispatch')
class OMIMAutocompleteView(AbstractOntologyTermAutocompleteView):
    def _get_ontology_service(self):
        return OntologyService.OMIM
