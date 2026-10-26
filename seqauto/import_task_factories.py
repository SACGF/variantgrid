"""
The DRAGEN TSO500 run-level files upload hands to seqauto: CombinedVariantOutput.tsv (the pair's analysis row) and
MetricsOutput.tsv (the run's per-library QC). Registered with upload through settings.IMPORT_TASK_FACTORY_IMPORTS.
"""
from django.conf import settings

from seqauto.models.models_upload import (
    UploadedDragenTSO500CombinedVariantOutput,
    UploadedDragenTSO500MetricsOutput,
)
from seqauto.tasks.import_dragen_tso500_combined_variant_output_task import (
    ImportDragenTSO500CombinedVariantOutputTask,
)
from seqauto.tasks.import_dragen_tso500_metrics_output_task import (
    ImportDragenTSO500MetricsOutputTask,
)
from seqauto.tso500 import dragen_metrics_output_parser
from upload.enums import UploadedFileTypes
from upload.import_task_factories.import_task_factory import ImportTaskFactory
from upload.tso500 import dragen_combined_variant_output_parser
from upload.upload_metadata import SEQUENCING_RUN


class DragenTSO500CombinedVariantOutputImportTaskFactory(ImportTaskFactory):
    """ Illumina DRAGEN TSO 500's CombinedVariantOutput.tsv - one vendor's format, not a standard.

        No variants come from it (the splice calls come from SpliceVariants.vcf, @see
        SpliceGirlImportTaskFactory): it is the pair's patient chain, seqauto links and TMB/MSI/GIS,
        which have no coordinates, so it is a single shot import rather than a VCF pipeline
        (@see seqauto.tasks.import_dragen_tso500_combined_variant_output_task).

        The file names its run 'NA', and that run keys the DragenTSO500CombinedVariantOutput row, so
        'sequencing_run' metadata is taken where the run cannot be found from the pair's sample IDs.

        A tsv full of gene symbols, so GeneListImportTaskFactory would otherwise claim it on its
        default ability of 1 """

    @property
    def enabled(self) -> bool:
        return settings.VARIANT_GENE_LEVEL_ENABLED

    def get_uploaded_file_type(self):
        return UploadedFileTypes.DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT

    def get_possible_extensions(self):
        return ['tsv']

    def get_metadata_keys(self):
        return frozenset({SEQUENCING_RUN})

    def get_data_classes(self):
        return [UploadedDragenTSO500CombinedVariantOutput]

    def get_processing_ability(self, user, filename, file_extension):
        if dragen_combined_variant_output_parser.can_process_file(filename):
            return 1000
        return 0

    def create_import_task(self, upload_pipeline):
        return ImportDragenTSO500CombinedVariantOutputTask.si(upload_pipeline.pk)


class DragenTSO500MetricsOutputImportTaskFactory(ImportTaskFactory):
    """ Illumina DRAGEN TSO 500's MetricsOutput.tsv - a run's per-library QC, one column per sample
        of every pair on the run. It has no variants and no coordinates, so it is a single shot
        import rather than a VCF pipeline (@see seqauto.tasks.import_dragen_tso500_metrics_output_task).

        The file names the run nowhere, and its column names alone do not identify a library across
        runs, so 'sequencing_run' is required metadata and is part of the LibraryQC key.

        A tsv, so GeneListImportTaskFactory would otherwise claim it on its default ability of 1 """

    @property
    def enabled(self) -> bool:
        return settings.VARIANT_GENE_LEVEL_ENABLED

    def get_uploaded_file_type(self):
        return UploadedFileTypes.DRAGEN_TSO500_METRICS_OUTPUT

    def get_possible_extensions(self):
        return ['tsv']

    def get_metadata_keys(self):
        return frozenset({SEQUENCING_RUN})

    def get_data_classes(self):
        return [UploadedDragenTSO500MetricsOutput]

    def get_processing_ability(self, user, filename, file_extension):
        if dragen_metrics_output_parser.can_process_file(filename):
            return 1000
        return 0

    def create_import_task(self, upload_pipeline):
        return ImportDragenTSO500MetricsOutputTask.si(upload_pipeline.pk)
