from typing import Optional

from django.contrib.auth.models import User
from guardian.shortcuts import assign_perm
import logging
from library.genomics.bed_file import BedFileReader
from snpdb.models import GenomicIntervalsCollection, ImportStatus, GenomicIntervalsCategory, GenomeBuild
from upload.models import UploadedBed
from upload.tasks.import_task import ImportTask
from variantgrid.celery import app


class ImportBedFileTask(ImportTask):
    """ It's not always possible to know what genome build a bed file is, so this may
        not be able to complete the processing (leaving at import_status REQUIRES_USER_INPUT) """

    def process_items(self, uploaded_file):
        logging.debug("ImportBedFileTask: process items")

        uploaded_category = GenomicIntervalsCategory.objects.get(name="Uploaded")
        user = uploaded_file.user
        genomic_intervals_collection = GenomicIntervalsCollection.objects.create(name=uploaded_file.name,
                                                                                 category=uploaded_category,
                                                                                 user=user,
                                                                                 import_status=ImportStatus.IMPORTING)
        assign_perm('snpdb.view_genomicintervalscollection', user, genomic_intervals_collection)
        uploaded_bed = UploadedBed.objects.create(uploaded_file=uploaded_file,
                                                  genomic_intervals_collection=genomic_intervals_collection)

        bed_file = uploaded_file.get_filename()
        genome_build = self._get_genome_build(user, bed_file)
        if genome_build:
            genomic_intervals_collection.genome_build = genome_build
            uploaded_bed.process_bed_file()
            processed_records = genomic_intervals_collection.processed_records
        else:
            genomic_intervals_collection.import_status = ImportStatus.REQUIRES_USER_INPUT
            processed_records = 0  # Ok, will load after user input
        genomic_intervals_collection.save()
        return processed_records

    @staticmethod
    def _get_genome_build(user: User, bed_filename) -> Optional[GenomeBuild]:
        """ Attempt to get from:
            * track (header)
            * filename
            * if user only has 1 genome build """
        bf_reader = BedFileReader(bed_filename)
        genome_build_name = bf_reader.get_genome_build_name()
        if genome_build_name:
            return GenomeBuild.get_name_or_alias(genome_build_name)

        genome_build = GenomeBuild.detect_from_filename(bed_filename)
        if genome_build:
            return genome_build

        # If server has exactly 1 genome build, use that
        genome_builds = list(GenomeBuild.builds_with_annotation())
        if len(genome_builds) == 1:
            return genome_builds[0]
        return None


ImportBedFileTask = app.register_task(ImportBedFileTask())  # @UndefinedVariable
