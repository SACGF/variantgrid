from typing import Any, Union, Type

from django.db.models import Q
from django.dispatch import receiver
from django.utils.safestring import SafeString

from genes.models import GeneSymbol
from ontology.models import OntologyTerm, OntologyService
from patients.models import Patient
from snpdb.search2 import SearchResponseRecordAbstract, search_signal, SearchInput, SearchResponse


class SearchResponsePatient(SearchResponseRecordAbstract[Patient]):

    @classmethod
    def result_class(cls) -> Type:
        return Patient


@receiver(search_signal, sender=SearchInput)
def patient_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    if search_input.matches_has_alpha():
        search_string = search_input.search_string
        parts = search_string.split(",")
        if len(parts) == 2:
            (last_name, first_name) = parts
            q_last = Q(last_name__iexact=last_name.strip())
            q_first = Q(first_name__iexact=first_name.strip())
            patient_q = q_last & q_first
        else:
            patient_q = Q(last_name__iexact=search_string) | Q(first_name__iexact=search_string)

        response: SearchResponse[SearchResponsePatient] = SearchResponse(SearchResponsePatient)
        response.extend(Patient.filter_for_user(search_input.user).filter(patient_q))
        return response
