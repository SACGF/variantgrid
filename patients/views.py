import mimetypes

import pandas as pd
from django.conf import settings
from django.http.response import HttpResponse, JsonResponse
from django.shortcuts import render, get_object_or_404
from django.views.decorators.http import require_POST
from jfu.http import upload_receive, UploadResponse, JFUResponse

from annotation.models.models_phenotype_match import TextPhenotypeMatch
from annotation.phenotype_matching import create_phenotype_description
from library.django_utils import add_save_message, set_form_read_only
from library.file_utils import rm_if_exists
from library.guardian_utils import assign_permission_to_user_and_groups
from library.log_utils import log_traceback
from library.utils import invert_dict
from ontology.forms import OMIMForm, HPOForm, HGNCForm, MONDOForm
from patients import forms
from patients.forms import PatientSearchForm, PatientContactForm
from patients.models import PatientColumns, PatientRecords, Patient, PatientModification, PatientRecordOriginType, \
    PatientAttachment
from snpdb.models import Sample
from uicore.utils.form_helpers import form_helper_horizontal


def view_patient(request, patient_id):
    patient = Patient.get_for_user(request.user, patient_id)
    form = forms.PatientForm(request.POST or None, instance=patient, user=request.user)
    form.helper = form_helper_horizontal()

    has_write_permission = patient.can_write(request.user)
    if not has_write_permission:
        set_form_read_only(form)
    show_read_only_patient_dob = not settings.PATIENTS_READ_ONLY_SHOW_AGE_NOT_DOB

    if request.method == "POST":
        valid = form.is_valid()
        if valid:
            patient = form.save()

        add_save_message(request, valid, "Patient")

    specimens = f" ({patient.num_specimens})" if patient.num_specimens else ''
    existing_files = get_patient_attachment_file_dicts(patient)

    context = {"patient": patient,
               "form": form,
               "specimens": specimens,
               "existing_files": existing_files,
               "show_read_only_patient_dob": show_read_only_patient_dob,
               "has_write_permission": has_write_permission}
    return render(request, 'patients/view_patient.html', context)


def view_patient_contact_tab(request, patient_id):
    patient = Patient.get_for_user(request.user, patient_id)
    form = PatientContactForm(request.POST or None, instance=patient)
    has_write_permission = patient.can_write(request.user)
    if not has_write_permission:
        set_form_read_only(form)

    if request.method == "POST":
        valid = form.is_valid()
        if valid:
            patient = form.save()
        add_save_message(request, valid, "Patient")

    context = {"patient": patient,
               "form": form,
               "has_write_permission": has_write_permission}
    return render(request, 'patients/view_patient_contact_tab.html', context)


def view_patient_specimens(request, patient_id):
    patient = Patient.get_for_user(request.user, patient_id)
    specimen_formset = forms.PatientSpecimenFormSet(request.POST or None, instance=patient)
    if request.method == "POST":
        valid = specimen_formset.is_valid()
        if valid:
            specimen_formset.save()
        add_save_message(request, valid, "Patient Specimen")

    context = {"patient": patient,
               "num_specimens": patient.num_specimens,
               "specimen_formset": specimen_formset,
               "has_write_permission": patient.can_write(request.user)}
    return render(request, 'patients/view_patient_specimens.html', context)


def view_patient_genes(request, patient_id):
    patient = Patient.get_for_user(request.user, patient_id)
    context = {"patient": patient}
    return render(request, 'patients/view_patient_genes.html', context)


def view_patient_modifications(request, patient_id):
    patient = Patient.get_for_user(request.user, patient_id)
    patient_modifications = patient.patientmodification_set.all().order_by("-date")
    context = {"patient_modifications": patient_modifications}
    return render(request, 'patients/view_patient_modifications.html', context)


@require_POST
def patient_file_upload(request, patient_id):
    try:
        patient = Patient.get_for_user(request.user, patient_id)
        patient.check_can_write(request.user)
        uploaded_file = upload_receive(request)
        patient_attachment = PatientAttachment.objects.create(patient=patient, file=uploaded_file)
        file_dict = patient_attachment.get_file_dict()
    except Exception as e:
        log_traceback()
        file_dict = {"error": str(e)}

    return UploadResponse(request, file_dict)


@require_POST
def patient_file_delete(request, pk):
    success = False
    try:
        patient_attachment = PatientAttachment.objects.get(pk=pk)
        patient_attachment.patient.check_can_write(request.user)
        rm_if_exists(patient_attachment.file.path)
        if patient_attachment.thumbnail_path:
            rm_if_exists(patient_attachment.thumbnail_path)
        patient_attachment.delete()
        success = True
    except PatientAttachment.DoesNotExist:
        pass

    return JFUResponse(request, success)


def get_patient_attachment_file_dicts(patient):
    file_dicts = []
    for patient_attachment in patient.patientattachment_set.all():
        file_dicts.append(patient_attachment.get_file_dict())

    file_dicts = list(reversed(file_dicts))  # JFU adds most recent at the end
    return file_dicts


def view_patient_file_attachment(request, patient_attachment_id, thumbnail=False):
    """ This is not done via static files, so we add security later """

    patient_attachment = get_object_or_404(PatientAttachment, pk=patient_attachment_id)
    patient_attachment.patient.check_can_view(request.user)

    if thumbnail and patient_attachment.thumbnail_path:
        filename = patient_attachment.thumbnail_path
    else:
        filename = patient_attachment.file.path

    content_type, _ = mimetypes.guess_type(filename)
    with open(filename, "rb") as f:
        image_data = f.read()
    return HttpResponse(image_data, content_type=content_type)


def view_patient_file_attachment_thumbnail(request, patient_attachment_id):
    return view_patient_file_attachment(request, patient_attachment_id, thumbnail=True)


def patient_record_imports(request):
    context = {}
    return render(request, 'patients/patient_record_imports.html', context)


def view_patient_records(request, patient_records_id):
    patient_records = get_object_or_404(PatientRecords, pk=patient_records_id)
    patient_records.uploaded_file.check_can_view(request.user)
    context = {"patient_records": patient_records}
    return render(request, 'patients/view_patient_records.html', context)


def import_patient_records_details(request):
    context = {"column_descriptions": PatientColumns.COLUMN_DETAILS}

    return render(request, 'patients/import_patient_records_details.html', context)


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

    empty_row = {c: '' for c in PatientColumns.COLUMNS}
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


def patients(request):
    form = forms.PatientForm(request.POST or None, user=request.user)
    if request.method == "POST":
        valid = form.is_valid()
        if valid:
            patient = form.save()

            assign_permission_to_user_and_groups(request.user, patient)
            form = forms.PatientForm(user=request.user)  # clear form for next patient
            msg = f"Patient #{patient.pk}: {str(patient)}"
        else:
            msg = "Patient"
        add_save_message(request, valid, msg, created=True)
        initially_hide_create_patient_form = False
    else:
        initially_hide_create_patient_form = True

    phenotype_match_graphs = TextPhenotypeMatch.objects.exists()
    context = {"phenotype_match_graphs": phenotype_match_graphs,
               "patient_search_form": PatientSearchForm(),
               'omim_form': OMIMForm(),
               'hpo_form': HPOForm(),
               'mondo_form': MONDOForm(),
               'hgnc_form': HGNCForm(),
               "initially_hide_create_patient_form": initially_hide_create_patient_form,
               "form": form}

    return render(request, 'patients/patients.html', context)


@require_POST
def phenotypes_matches(request):
    """ Live phenotype edits (unsaved changes) are sent here. """
    phenotype_text = request.POST['phenotype_text']
    phenotype_description = create_phenotype_description(phenotype_text)
    results = phenotype_description.get_results()
    phenotype_description.delete()
    return JsonResponse(results, safe=False)


def patient_term_matches(request):
    context = {}
    return render(request, 'patients/patient_term_matches.html', context)


def patient_term_approvals(request, patient_id_offset=0, num_patients_per_page=20, show_approved=False):
    filter_kwargs = {"patient_text_phenotype__phenotype_description__isnull": False}
    if not show_approved:
        filter_kwargs["patient_text_phenotype__approved_by__isnull"] = True

    if patient_id_offset:
        filter_kwargs["pk__gt"] = patient_id_offset
    patients_qs = Patient.objects.filter(**filter_kwargs).order_by("pk")
    num_total = patients_qs.count()

    end = min(num_patients_per_page, num_total)
    remaining_patients = num_total - end

    patients = list(patients_qs[:end])

    patient_results = {}
    max_patient_id = 0
    for p in patients:
        patient_results[p.pk] = p.patient_text_phenotype.phenotype_description.get_results()
        max_patient_id = p.pk

    context = {"patients": patients,
               "patient_results": patient_results,
               "remaining_patients": remaining_patients,
               "max_patient_id": max_patient_id,
               "show_approved": show_approved}
    return render(request, 'patients/patient_term_approvals.html', context)


def bulk_patient_term(request, patient_id_offset=0, num_patients_per_page=20):
    return patient_term_approvals(request, patient_id_offset=patient_id_offset, num_patients_per_page=num_patients_per_page, show_approved=True)


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
