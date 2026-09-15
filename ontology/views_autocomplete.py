import abc
import operator
from functools import reduce

from django.db.models import Case, IntegerField, Q, Value, When
from django.db.models.functions import Length
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page

from library.constants import HOUR_SECS
from library.django_utils.autocomplete_utils import AutocompleteView
from ontology.models import OntologyIdNormalized, OntologyService, OntologyTerm, OntologyTermStatus


class AbstractOntologyTermAutocompleteView(abc.ABC, AutocompleteView):
    fields = ['name']

    @abc.abstractmethod
    def _get_ontology_service(self):
        pass

    def sort_queryset(self, qs):
        return qs.order_by(Length("name").asc(), 'name')

    def _get_term_id_q(self) -> Q | None:
        """ A bare number matches the term with that index ("1061" -> HPO:0001061 in the HPO autocomplete),
            a prefixed id matches after normalisation ("HPO:1061", "mondo_7") """
        q = self.q.strip()
        if q.isdigit():
            return Q(index=int(q))
        try:
            return Q(id=OntologyIdNormalized.normalize(q).full_id)
        except ValueError:
            return None

    def get_queryset(self):
        user = self.request.user
        qs = self.get_user_queryset(user)
        if not user.is_authenticated:
            return qs.none()

        if self.q:
            name_q = Q(name__icontains=self.q)
            if term_id_q := self._get_term_id_q():
                # The term whose id was typed goes first, ahead of any names that happen to contain the digits
                is_term_match = Case(When(term_id_q, then=Value(0)), default=Value(1), output_field=IntegerField())
                qs = qs.filter(name_q | term_id_q).annotate(is_term_match=is_term_match)
                return qs.order_by("is_term_match", Length("name").asc(), 'name')
            qs = qs.filter(name_q)

        return self.sort_queryset(qs)

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
        value = self.forwarded.get('ontology_service')
        if value is None:
            return None
        try:
            return OntologyService(value)
        except ValueError:
            return None


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
