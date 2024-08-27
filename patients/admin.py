from django.contrib import admin

from patients import models
from snpdb.admin_utils import ModelAdminBasics


@admin.register(models.Patient)
class PatientAdmin(ModelAdminBasics):
    pass


@admin.register(models.Tissue)
class TissueAdmin(ModelAdminBasics):
    pass


@admin.register(models.Specimen)
class SpecimenAdmin(ModelAdminBasics):
    pass


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
