from django.contrib import admin

from patients import models
from snpdb.admin_utils import ModelAdminBasics


@admin.register(models.Patient)
class PatientAdmin(ModelAdminBasics):
    list_display = ("id", "patient_code", "family_code", "first_name", "last_name", "sex", "modified")
    search_fields = ("patient_code", "family_code", "first_name", "last_name")


@admin.register(models.Tissue)
class TissueAdmin(ModelAdminBasics):
    pass


@admin.register(models.Specimen)
class SpecimenAdmin(ModelAdminBasics):
    list_display = ("id", "reference_id", "patient", "tissue", "mutation_type", "collection_date")
    search_fields = ("reference_id",)


@admin.register(models.Extraction)
class ExtractionAdmin(ModelAdminBasics):
    list_display = ("id", "specimen", "reference_id", "nucleic_acid_source")
    search_fields = ("reference_id", "specimen__reference_id")


@admin.register(models.PatientImport)
class PatientImportAdmin(ModelAdminBasics):
    pass


@admin.register(models.PatientModification)
class PatientModificationAdmin(ModelAdminBasics):
    pass


@admin.register(models.PatientComment)
class PatientCommentAdmin(ModelAdminBasics):
    pass


@admin.register(models.Clinician)
class ClinicianAdmin(ModelAdminBasics):
    pass


@admin.register(models.PatientRecords)
class PatientRecordsAdmin(ModelAdminBasics):
    pass


@admin.register(models.PatientRecord)
class PatientRecordAdmin(ModelAdminBasics):
    pass
