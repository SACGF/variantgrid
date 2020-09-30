import logging

from patients.import_records import import_patient_records
from patients.models import PatientRecords, PatientImport
from upload.models import UploadedPatientRecords
from upload.tasks.import_task import ImportTask
from variantgrid.celery import app


class ImportPatientRecords(ImportTask):
    def process_items(self, uploaded_file):
        patient_import = PatientImport.objects.create()
        patient_records = PatientRecords.objects.create(patient_import=patient_import)
        logging.info("Created patient_records: %s", patient_records)
        UploadedPatientRecords.objects.create(uploaded_file=uploaded_file,
                                              patient_records=patient_records)

        return import_patient_records(patient_records)

ImportPatientRecords = app.register_task(ImportPatientRecords())  # @UndefinedVariable
