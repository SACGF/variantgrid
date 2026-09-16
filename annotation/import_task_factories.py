import cyvcf2

from annotation.tasks.import_clinvar_vcf_task import (
    ImportClinVarSuccessTask,
    ImportCreateVersionForClinVarVCFTask,
    ProcessClinVarVCFDataTask,
)
from annotation.vcf_files.import_clinvar_vcf import check_can_import_clinvar
from library.genomics.vcf_utils import cyvcf2_header_get
from upload.import_task_factories.abstract_vcf_import_task_factory import (
    AbstractVCFImportTaskFactory,
)
from upload.models import UploadedClinVarVersion, UploadedFileTypes


class ImportClinVarTaskFactory(AbstractVCFImportTaskFactory):
    def get_uploaded_file_type(self):
        return UploadedFileTypes.CLINVAR

    def get_possible_extensions(self):
        return ['vcf']

    def get_data_classes(self):
        return [UploadedClinVarVersion]

    def get_processing_ability(self, user, filename, file_extension, **kwargs):
        # cyvcf2 is what the pipeline reads with, and unlike PyVCF it tolerates the extra keys some callers
        # put in FILTER lines (eg Illumina's Number=1,Type=String) - this runs for every uploaded VCF
        reader = cyvcf2.VCF(filename)
        if cyvcf2_header_get(reader, "source") == "ClinVar":
            check_can_import_clinvar(user)
            return 2  # more specific than normal VCF
        return 0

    def get_create_data_from_vcf_header_task_class(self):
        return ImportCreateVersionForClinVarVCFTask

    def get_known_variants_parallel_vcf_processing_task_class(self):
        return ProcessClinVarVCFDataTask

    def get_finish_task_classes(self):
        return [ImportClinVarSuccessTask]
