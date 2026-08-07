import logging

from django.conf import settings

from pedigree.ped.import_ped import automatch_pedigree_samples, import_ped
from upload.models import UploadedPedFile
from upload.tasks.import_task import ImportTask
from variantgrid.celery import app


class ImportPedTask(ImportTask):

    def process_items(self, file_upload):
        logging.debug("ImportPedTask: process items")

        ped_file, families = import_ped(file_upload.get_file(), file_upload.name, file_upload.user)

        min_matching_samples = getattr(settings, "PEDIGREE_MIN_COHORT_SAMPLE_MATCHES_FOR_AUTO_MATCH", None)
        if min_matching_samples:
            automatch_pedigree_samples(file_upload.user, families, min_matching_samples)

        UploadedPedFile.objects.create(file_upload=file_upload,
                                       ped_file=ped_file)

        inserted = sum([p.get_records_count() for p in families])
        return inserted


ImportPedTask = app.register_task(ImportPedTask())  # @UndefinedVariable
