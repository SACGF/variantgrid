from django.db.models import Q
from patients.models import Patient
from snpdb.search2 import search_receiver, SearchInputInstance, SearchExample


@search_receiver(
    search_type=Patient,
    example=SearchExample(
        note="(Last Name, First Name) or Last Name or First Name",
        example="Smith, Alvin"
    )
)
def patient_search(search_input: SearchInputInstance):
    search_string = search_input.search_string
    parts = search_string.split(",")
    if len(parts) == 2:
        (last_name, first_name) = parts
        q_last = Q(last_name__iexact=last_name.strip())
        q_first = Q(first_name__iexact=first_name.strip())
        patient_q = q_last & q_first
    else:
        patient_q = Q(last_name__iexact=search_string) | Q(first_name__iexact=search_string)

    yield Patient.filter_for_user(search_input.user).filter(patient_q)
