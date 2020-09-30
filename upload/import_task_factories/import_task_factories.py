from django.conf import settings

from library.utils import import_class, full_class_name
from patients.models import PatientColumns
from upload.import_task_factories.abstract_vcf_import_task_factory import AbstractVCFImportTaskFactory
from upload.import_task_factories.import_task_factory import ImportTaskFactory
from upload.models import UploadedFileTypes, UploadedBed, UploadedExpressionFile, \
    UploadedGeneList, UploadedPatientRecords, UploadedPedFile, UploadedVCF, \
    UploadedGeneCoverage, UploadStep
from upload.models_enums import UploadStepTaskType
from upload.tasks.import_bedfile_task import ImportBedFileTask
from upload.tasks.import_expression_task import ImportExpressionTask
from upload.tasks.import_gene_coverage_task import ImportGeneCoverageTask
from upload.tasks.import_gene_list_task import ImportGeneListTask
from upload.tasks.import_patient_records_task import ImportPatientRecords
from upload.tasks.import_ped_task import ImportPedTask
from upload.tasks.vcf.genotype_vcf_tasks import VCFCheckAnnotationTask, \
    ProcessGenotypeVCFDataTask, ImportGenotypeVCFSuccessTask, \
    UpdateVariantZygosityCountsTask, SampleLocusCountsTask, \
    ImportCreateVCFModelForGenotypeVCFTask
from upload.tasks.vcf.import_vcf_tasks import ProcessVCFSetMaxVariantTask, \
    ImportCreateUploadedVCFTask, ProcessVCFLinkAllelesSetMaxVariantTask, LiftoverCompleteTask, LiftoverCreateVCFTask
import pandas as pd


class BedImportTaskFactory(ImportTaskFactory):

    def get_uploaded_file_type(self):
        return UploadedFileTypes.BED

    def get_possible_extensions(self):
        return ['bed']

    def get_data_classes(self):
        return [UploadedBed]

    def create_import_task(self, upload_pipeline):
        return ImportBedFileTask.si(upload_pipeline.pk)


class ExpressionImportTaskFactory(ImportTaskFactory):

    def get_uploaded_file_type(self):
        return UploadedFileTypes.CUFFDIFF

    def get_possible_extensions(self):
        return ['diff']

    def get_data_classes(self):
        return [UploadedExpressionFile]

    def create_import_task(self, upload_pipeline):
        return ImportExpressionTask.si(upload_pipeline.pk)


class GeneListImportTaskFactory(ImportTaskFactory):

    def get_uploaded_file_type(self):
        return UploadedFileTypes.GENE_LIST

    def get_possible_extensions(self):
        return ['tsv', 'csv', 'txt']

    def get_data_classes(self):
        return [UploadedGeneList]

    def get_processing_ability(self, user, filename, file_extension):
        # TODO: Check actually contains genes... but will be overriden by other CSV handlers
        return 1

    def create_import_task(self, upload_pipeline):
        return ImportGeneListTask.si(upload_pipeline.pk)


class GeneCoverageImportTaskFactory(ImportTaskFactory):

    def get_uploaded_file_type(self):
        return UploadedFileTypes.GENE_COVERAGE

    def get_possible_extensions(self):
        return ['txt']

    def get_data_classes(self):
        return [UploadedGeneCoverage]

    def get_processing_ability(self, user, filename, file_extension):
        if ImportGeneCoverageTask.can_process_file(filename):
            return 1000000
        return 0

    def create_import_task(self, upload_pipeline):
        return ImportGeneCoverageTask.si(upload_pipeline.pk)


class PatientRecordsImportTaskFactory(ImportTaskFactory):

    def get_uploaded_file_type(self):
        return UploadedFileTypes.PATIENT_RECORDS

    def get_possible_extensions(self):
        return ['csv']

    def get_data_classes(self):
        return [UploadedPatientRecords]

    def get_processing_ability(self, user, filename, file_extension):
        df = pd.read_csv(filename, encoding='unicode_escape')
        if PatientColumns.PATIENT_LAST_NAME in df.columns:
            return 1000

    def create_import_task(self, upload_pipeline):
        return ImportPatientRecords.si(upload_pipeline.pk)


class PedImportTaskFactory(ImportTaskFactory):

    def get_uploaded_file_type(self):
        return UploadedFileTypes.PED

    def get_possible_extensions(self):
        return ['ped']

    def get_data_classes(self):
        return [UploadedPedFile]

    def create_import_task(self, upload_pipeline):
        return ImportPedTask.si(upload_pipeline.pk)


class GenotypeVCFImportFactory(AbstractVCFImportTaskFactory):
    """ This doesn't do anything, just needed to be consistent with the other ones... """

    def get_uploaded_file_type(self):
        return UploadedFileTypes.VCF

    def get_data_classes(self):
        return [UploadedVCF]

    def get_create_data_from_vcf_header_task_class(self):
        return ImportCreateVCFModelForGenotypeVCFTask

    def get_known_variants_parallel_vcf_processing_task_class(self):
        return ProcessGenotypeVCFDataTask

    def get_post_data_insertion_classes(self):
        return [VCFCheckAnnotationTask,
                UpdateVariantZygosityCountsTask,
                SampleLocusCountsTask]

    def get_finish_task_classes(self):
        finish_task_classes = []

        for class_name in settings.FINISH_IMPORT_VCF_STEP_TASKS_CLASSES:
            finish_task_classes.append(import_class(class_name))

        finish_task_classes.append(ImportGenotypeVCFSuccessTask)
        return finish_task_classes


class VCFInsertVariantsOnlyImportFactory(AbstractVCFImportTaskFactory):

    def get_uploaded_file_type(self):
        return UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY

    def get_data_classes(self):
        return [UploadedVCF]

    def get_processing_ability(self, user, filename, file_extension):
        # Variants Only should only ever be explicitly chosen as a file type
        # Never detected and set on an uploaded file
        return 0

    def get_create_data_from_vcf_header_task_class(self):
        return ImportCreateUploadedVCFTask

    def get_known_variants_parallel_vcf_processing_task_class(self):
        return ProcessVCFSetMaxVariantTask

    def get_post_data_insertion_classes(self):
        return [VCFCheckAnnotationTask]


class LiftoverImportFactory(AbstractVCFImportTaskFactory):
    """ Converts variants to a different genome build, linking via Allele records
        The Allele PK is written into VCF record column #3 (ID) which BulkAlleleLinkingVCFProcessor uses to
        links the new variant to Allele via VariantAllele @see https://github.com/SACGF/variantgrid/wiki/Liftover """
    def get_uploaded_file_type(self):
        return UploadedFileTypes.LIFTOVER

    def get_data_classes(self):
        return [UploadedVCF]

    def get_processing_ability(self, user, filename, file_extension):
        # Variants Only should only ever be explicitly chosen as a file type
        # Never detected and set on an uploaded file
        return 0

    def get_pre_vcf_task(self, upload_pipeline):
        # If source_vcf/source_genome_build are set, need to convert to produce allele VCF for genome_build
        liftover = upload_pipeline.uploaded_file.uploadedliftover.liftover
        pre_vcf_task = None
        if liftover.source_genome_build:
            output_filename = upload_pipeline.uploaded_file.get_filename()
            unknown_variants_step = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                                              name="Liftover Source VCF",
                                                              sort_order=self.get_sort_order(),
                                                              task_type=UploadStepTaskType.CELERY,
                                                              script=full_class_name(LiftoverCreateVCFTask),
                                                              input_filename=liftover.source_vcf,
                                                              output_filename=output_filename)
            pre_vcf_task = LiftoverCreateVCFTask.si(unknown_variants_step.pk, 0)
        return pre_vcf_task

    def get_create_data_from_vcf_header_task_class(self):
        return ImportCreateUploadedVCFTask

    def get_known_variants_parallel_vcf_processing_task_class(self):
        return ProcessVCFLinkAllelesSetMaxVariantTask

    def get_post_data_insertion_classes(self):
        return [VCFCheckAnnotationTask]

    def get_finish_task_classes(self):
        task_classes = super().get_finish_task_classes()
        return [LiftoverCompleteTask] + task_classes
