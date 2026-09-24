from django.contrib.auth.models import User
from django.db.models import CharField, Count, Value
from django.db.models.functions import Concat, Lower
from django.dispatch import receiver

from annotation.models.models_phenotype_match import (
    PHENOTYPE_ONTOLOGY_SERVICE_LABELS,
    patient_phenotype_terms,
    patients_qs_for_ontology_term,
)
from library.preview_request import PreviewKeyValue, preview_extra_signal
from ontology.models import OntologyTerm
from patients.models import Patient
from snpdb.models import Sample
from snpdb.search import HAS_3_ANY, SearchExample, SearchInputInstance, search_receiver


@search_receiver(
    search_type=Patient,
    pattern=HAS_3_ANY,
    example=SearchExample(
        note="3 or more characters in patient's name or patient code",
        examples=["Smith Alvin"]
    )
)
def patient_search(search_input: SearchInputInstance):
    qs = Patient.objects.annotate(
        combined_name=Concat('first_name', Value(' '), 'last_name', output_field=CharField())
    )
    name_q = search_input.q_words('combined_name')
    code_q = search_input.q_words('patient_code')
    yield qs.filter(name_q | code_q).order_by(Lower('last_name'), Lower('first_name'))


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


@receiver(preview_extra_signal, sender=Sample)
def sample_preview_patient_extra(sender, user: User, obj: Sample, **kwargs):
    """ The sample's patient and the phenotype terms matched in their record, one row per ontology.
        Only for a patient the user can view - as the sample gene matrix does """
    if not obj.patient_id:
        return None
    if not (patient := Patient.filter_for_user(user).filter(pk=obj.patient_id).first()):
        return None
    patient_row = PreviewKeyValue(key="Patient", value=patient.display_identity, icon=Patient.preview_icon())
    return [patient_row, *_phenotype_preview_rows(patient)]


@receiver(preview_extra_signal, sender=Patient)
def patient_preview_phenotype_extra(sender, user: User, obj: Patient, **kwargs):
    return _phenotype_preview_rows(obj)


def _phenotype_preview_rows(patient: Patient) -> list[PreviewKeyValue]:
    """ The terms matched in the patient's phenotype text, one row per ontology """
    rows = []
    if phenotype_terms := patient_phenotype_terms([patient]).get(patient.pk):
        for service_label in PHENOTYPE_ONTOLOGY_SERVICE_LABELS.values():
            if terms := phenotype_terms.terms.get(service_label):
                rows.append(PreviewKeyValue(key=service_label, value=", ".join(term.name for term in terms)))
    return rows
