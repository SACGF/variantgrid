import vcf

from annotation.tasks.import_clinvar_vcf_task import ImportCreateVersionForClinVarVCFTask, ProcessClinVarVCFDataTask, \
    ImportClinVarSuccessTask
from annotation.vcf_files.import_clinvar_vcf import check_can_import_clinvar
from upload.import_task_factories.abstract_vcf_import_task_factory import AbstractVCFImportTaskFactory
from upload.models import UploadedFileTypes


class ImportClinVarTaskFactory(AbstractVCFImportTaskFactory):
    def get_uploaded_file_type(self):
        return UploadedFileTypes.CLINVAR

    def get_possible_extensions(self):
        return ['vcf']

    def get_data_classes(self):
        return []  # TODO?

    def get_processing_ability(self, user, filename, file_extension, **kwargs):
        # Even if it originally was a .gz, will be using decompressed stream here
        reader = vcf.Reader(filename=filename)
        source = reader.metadata.get("source")
        if source and source == ['ClinVar']:
            check_can_import_clinvar(user)
            return 2  # more specific than normal VCF
        return 0

    def get_create_data_from_vcf_header_task_class(self):
        return ImportCreateVersionForClinVarVCFTask

    def get_known_variants_parallel_vcf_processing_task_class(self):
        return ProcessClinVarVCFDataTask

    def get_finish_task_classes(self):
        return [ImportClinVarSuccessTask]
