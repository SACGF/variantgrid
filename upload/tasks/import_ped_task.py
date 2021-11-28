import logging

from django.conf import settings

from pedigree.ped.import_ped import import_ped, automatch_pedigree_samples
from upload.models import UploadedPedFile
from upload.tasks.import_task import ImportTask
from variantgrid.celery import app


class ImportPedTask(ImportTask):

    def process_items(self, uploaded_file):
        logging.debug("ImportPedTask: process items")

        ped_file, families = import_ped(uploaded_file.get_file(), uploaded_file.name, uploaded_file.user)

        min_matching_samples = getattr(settings, "PEDIGREE_MIN_COHORT_SAMPLE_MATCHES_FOR_AUTO_MATCH", None)
        if min_matching_samples:
            automatch_pedigree_samples(uploaded_file.user, families, min_matching_samples)

        UploadedPedFile.objects.create(uploaded_file=uploaded_file,
                                       ped_file=ped_file)

        inserted = sum([p.get_records_count() for p in families])
        return inserted


ImportPedTask = app.register_task(ImportPedTask())  # @UndefinedVariable
