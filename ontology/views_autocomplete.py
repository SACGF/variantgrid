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

    def _get_term_id_q(self) -> tuple[Q, Q] | None:
        """ Returns (id_q, exact_q): digits match anywhere in the id ("123" -> HP:0000123, HP:0012323),
            a prefixed id ("HPO:123", "mondo_7") the same within that prefix. exact_q is the term whose
            index is exactly the digits typed, for ranking """
        q = self.q.strip()
        if q.isdigit():
            return Q(id__contains=q), Q(index=int(q))
        try:
            normalized = OntologyIdNormalized.normalize(q)
        except ValueError:
            return None
        digits = normalized.postfix.lstrip("0") or "0"
        id_q = Q(id__startswith=f"{normalized.prefix}:") & Q(id__icontains=digits)
        return id_q, Q(id=normalized.full_id)

    def get_queryset(self):
        user = self.request.user
        qs = self.get_user_queryset(user)
        if not user.is_authenticated:
            return qs.none()

        if words := self.q.split():
            # Every word, in any order: "muscular dystrophy duchenne" finds "Duchenne muscular dystrophy"
            name_q = reduce(operator.and_, (Q(name__icontains=word) for word in words))
            if term_id_qs := self._get_term_id_q():
                id_q, exact_q = term_id_qs
                # The term with exactly that id first, then other ids containing the digits, then names that do
                term_match_rank = Case(When(exact_q, then=Value(0)), When(id_q, then=Value(1)), default=Value(2),
                                       output_field=IntegerField())
                qs = qs.filter(name_q | id_q).annotate(term_match_rank=term_match_rank)
                return qs.order_by("term_match_rank", Length("name").asc(), 'name')
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
