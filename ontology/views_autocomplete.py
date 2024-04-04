import abc
import operator
from functools import reduce

from django.db.models import Q
from django.db.models.functions import Length
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page

from library.constants import HOUR_SECS
from library.django_utils.autocomplete_utils import AutocompleteView
from ontology.models import OntologyTerm, OntologyService, OntologyTermStatus


class AbstractOntologyTermAutocompleteView(abc.ABC, AutocompleteView):
    fields = ['name']

    @abc.abstractmethod
    def _get_ontology_service(self):
        pass

    def sort_queryset(self, qs):
        return qs.order_by(Length("name").asc(), 'name')

    def get_user_queryset(self, user):
        qs = OntologyTerm.objects.all()
        filters = []
        if ontology_service := self._get_ontology_service():
            if ontology_service in [OntologyService.HPO, OntologyService.MONDO, OntologyService.OMIM]:
                # These we want to hide obsolete/gene etc
                filters.extend([
                    Q(status=OntologyTermStatus.CONDITION),
                    Q(name__isnull=False),
                    ~Q(name=''),
                ])
            filters.append(Q(ontology_service=ontology_service))

        if filters:
            q = reduce(operator.and_, filters)
            qs = qs.filter(q)

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


@method_decorator(cache_page(HOUR_SECS), name='dispatch')
class HGNCAutocompleteView(AbstractOntologyTermAutocompleteView):
    def _get_ontology_service(self):
        return OntologyService.HGNC


@method_decorator(cache_page(HOUR_SECS), name='dispatch')
class MONDOAutocompleteView(AbstractOntologyTermAutocompleteView):
    def _get_ontology_service(self):
        return OntologyService.MONDO

    def get_user_queryset(self, user):
        qs = super().get_user_queryset(user)
        if self.forwarded.get('gene_disease'):
            qs = qs.filter(subject__extra__strongest_classification__isnull=False).distinct()
        return qs
