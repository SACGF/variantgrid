from django.contrib import admin

from patients import models

admin.site.register(models.Patient)
admin.site.register(models.Tissue)
admin.site.register(models.Specimen)
admin.site.register(models.PatientImport)
admin.site.register(models.PatientModification)
admin.site.register(models.PatientComment)
admin.site.register(models.Clinician)
admin.site.register(models.PatientRecords)
admin.site.register(models.PatientRecord)
