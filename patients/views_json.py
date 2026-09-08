import pandas as pd
from django.http import Http404
from django.http.response import HttpResponse, JsonResponse
from django.views.decorators.http import require_POST

from annotation.phenotype_matching import create_phenotype_description
from library.utils import invert_dict
from patients import forms
from patients.models import (
    Patient,
    PatientColumns,
    PatientModification,
    PatientRecordOriginType,
)
from patients.sample_grouping import SOURCE_LEVELS, get_patient_sample_tree
from snpdb.models import GenomeBuild, Sample


def _get_sample_source(user, level: str, pk: int):
    """ A level of Patient -> Specimen -> Extraction -> Sample, loaded through its own permissions """
    source_level = SOURCE_LEVELS.get(level)
    if source_level is None:
        raise Http404(f"Unknown sample source level '{level}'")
    return source_level.model.get_for_user(user, pk)


def _get_request_genome_build(request):
    if genome_build_name := request.GET.get("genome_build"):
        return GenomeBuild.get_name_or_alias(genome_build_name)
    return None


def sample_group_tree(request, level, pk):
    """ The whole patient the object belongs to, with the picked row's subtree flagged - what the
        SampleNode editor draws so moving up or down a level doesn't need another search """
    source = _get_sample_source(request.user, level, pk)
    tree = get_patient_sample_tree(request.user, level, source, _get_request_genome_build(request))
    return JsonResponse(tree)


def example_upload_csv_empty(request):
    """ headers only """
    sample_qs = Sample.objects.none()
    filename = "example_patient_upload"
    return get_patient_upload_csv(filename, sample_qs)


def example_upload_csv_all(request):
    sample_qs = Sample.filter_for_user(request.user)
    columns_lookup = invert_dict(PatientColumns.SAMPLE_QUERYSET_PATH)
    filename = f"{request.user}_all_samples_upload"
    return get_patient_upload_csv(filename, sample_qs, columns_lookup=columns_lookup)


def example_upload_csv_no_patients(request):
    sample_qs = Sample.filter_for_user(request.user).filter(patient__isnull=True)
    filename = f"{request.user}_samples_without_patients_upload"
    return get_patient_upload_csv(filename, sample_qs)


def get_patient_upload_csv(filename, sample_qs, columns_lookup=None):
    response = HttpResponse(content_type='text/csv')
    filename = f"{filename}.csv"
    response['Content-Disposition'] = f'attachment; filename="{filename}"'

    if columns_lookup is None:
        columns_lookup = {
            "pk": PatientColumns.SAMPLE_ID,
            "name": PatientColumns.SAMPLE_NAME,
        }
    sample_values_qs = sample_qs.values(*columns_lookup)

    empty_row = dict.fromkeys(PatientColumns.COLUMNS, '')
    rows = []
    for values in sample_values_qs:
        data = empty_row.copy()
        for from_col, to_col in columns_lookup.items():
            val = values.get(from_col)
            if val is not None:
                data[to_col] = val
        rows.append(data)

    if not rows:
        rows = [empty_row]

    df = pd.DataFrame.from_dict(rows)
    df = df[PatientColumns.COLUMNS]
    response.write(df.to_csv(index=False))
    return response


@require_POST
def create_patient(request):
    form = forms.PatientForm(request.POST, user=request.user)
    valid = form.is_valid()
    data = {}
    if valid:
        patient = form.save()
        data["patient_id"] = patient.pk
        data["__str__"] = str(patient)
    else:
        data["error"] = form.errors

    return JsonResponse(data)


@require_POST
def phenotypes_matches(request):
    """ Live phenotype edits (unsaved changes) are sent here. """
    phenotype_text = request.POST['phenotype_text']
    phenotype_description = create_phenotype_description(phenotype_text)
    results = phenotype_description.get_results()
    phenotype_description.delete()
    return JsonResponse(results, safe=False)


@require_POST
def approve_patient_term(request):
    patient_id = request.POST["patient_id"]
    patient = Patient.get_for_user(request.user, patient_id)
    patient.patient_text_phenotype.approved_by = request.user
    patient.patient_text_phenotype.save()

    PatientModification.objects.create(patient=patient,
                                       user=request.user,
                                       description="Approved phenotype text.",
                                       origin=PatientRecordOriginType.MANUAL_VG_GUI)
    return HttpResponse()
