import logging

from annotation.models.models import ClinVarVersion
from annotation.vcf_files.import_clinvar_vcf import import_clinvar_vcf, \
    check_can_import_clinvar
from upload.models import UploadedClinVarVersion
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.vcf.vcf_import import vcf_detect_genome_build
from variantgrid.celery import app


class ImportCreateVersionForClinVarVCFTask(ImportVCFStepTask):
    """ Create Data from header """

    def process_items(self, upload_step):
        check_can_import_clinvar(upload_step.file_upload.user)
        filename = upload_step.input_filename
        genome_build = vcf_detect_genome_build(filename)

        upload_step.file_upload.store_sha256_hash()
        kwargs = {"sha256_hash": upload_step.file_upload.sha256_hash,
                  "genome_build": genome_build}
        clinvar_version, created = ClinVarVersion.objects.get_or_create(**kwargs)
        if created:
            logging.info("Created clinvar: %s", clinvar_version)
            clinvar_version.filename = filename
            clinvar_version.annotation_date = ClinVarVersion.get_annotation_date_from_filename(filename)
            clinvar_version.save()
        else:
            logging.info("Version already exists, deleting previous annotation")
            clinvar_version.delete_related_objects()
            clinvar_version.create_partition()  # Put it back...

        UploadedClinVarVersion.objects.get_or_create(file_upload=upload_step.file_upload,
                                                     clinvar_version=clinvar_version)
        return 0


class ProcessClinVarVCFDataTask(ImportVCFStepTask):
    """ This is called after any unknown variants have been inserted """

    def process_items(self, upload_step):
        return import_clinvar_vcf(upload_step)


class ImportClinVarSuccessTask(ImportVCFStepTask):

    def process_items(self, upload_step):
        # Cohort stats referencing the old ClinVar version are recomputed lazily when next read
        return 0


ImportCreateVersionForClinVarVCFTask = app.register_task(ImportCreateVersionForClinVarVCFTask())
ProcessClinVarVCFDataTask = app.register_task(ProcessClinVarVCFDataTask())
ImportClinVarSuccessTask = app.register_task(ImportClinVarSuccessTask())
