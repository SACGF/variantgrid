from patients.import_records import import_patient_records
from patients.models import PatientImport, PatientRecords
from upload.models import UploadedPatientRecords
from upload.tasks.import_task import ImportTask
from variantgrid.celery import app


class ImportPatientRecords(ImportTask):
    def process_items(self, file_upload):
        patient_import = PatientImport.objects.create()
        patient_records = PatientRecords.objects.create(patient_import=patient_import)
        UploadedPatientRecords.objects.create(file_upload=file_upload,
                                              patient_records=patient_records)

        return import_patient_records(patient_records)


ImportPatientRecords = app.register_task(ImportPatientRecords())  # @UndefinedVariable
