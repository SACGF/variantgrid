from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from upload.models import UploadedFileTypes
from upload.tasks.import_patient_records_task import ImportPatientRecords
from upload.upload_processing import create_upload_pipeline


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--user', required=True)
        parser.add_argument('patient_csv')

    def handle(self, *args, **options):
        username = options["user"]
        filename = options["patient_csv"]
        user = User.objects.get(username=username)
        ufpj = create_upload_pipeline(user, filename, UploadedFileTypes.GENE_LIST)

        task = ImportPatientRecords.si(ufpj.pk)
        task.apply()
