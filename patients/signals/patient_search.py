from django.contrib.auth.models import User
from django.db.models import Value, CharField, Count
from django.db.models.functions import Concat, Lower
from django.dispatch import receiver

from annotation.models import patients_qs_for_ontology_term
from library.preview_request import preview_extra_signal, PreviewKeyValue
from ontology.models import OntologyTerm
from patients.models import Patient
from snpdb.models import Sample
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, HAS_3_ANY


@search_receiver(
    search_type=Patient,
    pattern=HAS_3_ANY,
    example=SearchExample(
        note="3 or more characters in patient's name",
        examples=["Smith Alvin"]
    )
)
def patient_search(search_input: SearchInputInstance):
    yield Patient.objects.annotate(
        combined_name=Concat('first_name', Value(' '), 'last_name', output_field=CharField())
    ).filter(search_input.q_words('combined_name')).order_by(Lower('last_name'), Lower('first_name'))


@receiver(preview_extra_signal, sender=OntologyTerm)
def ontology_preview_patient_sample_extra(sender, user: User, obj: OntologyTerm, **kwargs):
    if not Patient.preview_enabled():
        return
    patients_qs = patients_qs_for_ontology_term(user, obj)
    data = patients_qs.aggregate(num_patients=Count("id", distinct=True), num_samples=Count("sample", distinct=True))
    extras = []
    if num_patients := data.get("num_patients"):
        extras.append(PreviewKeyValue.count(Patient, num_patients))
        if num_samples := data.get("num_samples"):
            extras.append(PreviewKeyValue.count(Sample, num_samples))
    return extras
