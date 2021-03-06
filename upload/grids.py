from django.shortcuts import get_object_or_404

from annotation.annotation_version_querysets import get_queryset_for_latest_annotation_version
from library.jqgrid_user_row_config import JqGridUserRowConfig
from snpdb.models.models_variant import Variant
from upload.models import UploadStep, ModifiedImportedVariant, UploadPipeline


class UploadStepsGrid(JqGridUserRowConfig):
    model = UploadStep
    caption = 'Upload Steps'
    fields = ['sort_order', 'id', 'name', 'status', 'items_to_process', 'items_processed', 'error_message', 'input_filename', 'output_filename', 'start_date', 'end_date', 'task_type', 'pipeline_stage', 'pipeline_stage_dependency', 'script', 'output_text', 'tool_version__name', 'tool_version__version', 'import_variant_table', 'celery_task']
    colmodel_overrides = {'name': {'width': 300},
                          'tool__name': {'label': 'Tool'},
                          'tool__version': {'label': 'Tool Version'}}

    def __init__(self, user, upload_pipeline_id):
        super().__init__(user)
        upload_pipeline = get_object_or_404(UploadPipeline, pk=upload_pipeline_id)
        upload_pipeline.uploaded_file.check_can_view(user)
        self.queryset = self.get_queryset(None).filter(upload_pipeline=upload_pipeline)  # .exclude(status=ProcessingStatus.SKIPPED)
        self.extra_config.update({'sortname': 'sort_order',
                                  'sortorder': 'desc'})


class UploadPipelineSkippedAnnotationGrid(JqGridUserRowConfig):
    model = Variant
    caption = 'Skipped Annotation'
    fields = ["id", "variantannotation__vep_skipped_reason", "variantannotation__annotation_run"]

    colmodel_overrides = {"id": {"hidden": True},
                          "variantannotation__annotation_run": {'formatter': 'formatAnnotationRunLink'}}

    def __init__(self, user, upload_pipeline_id):
        super().__init__(user)
        upload_pipeline = get_object_or_404(UploadPipeline, pk=upload_pipeline_id)
        vcf = upload_pipeline.uploadedvcf.vcf

        qs = get_queryset_for_latest_annotation_version(self.model, upload_pipeline.genome_build)
        qs = vcf.get_variant_qs(qs).filter(variantannotation__vep_skipped_reason__isnull=False)
        qs = Variant.annotate_variant_string(qs)

        field_names = list(self.get_field_names())
        field_names.insert(1, "variant_string")

        self.queryset = qs.values(*field_names)
        self.extra_config.update({'sortname': 'variant_string',
                                  'sortorder': 'asc'})

    def get_colmodels(self, remove_server_side_only=False):
        before_colmodels = [{'index': 'variant_string', 'name': 'variant_string', 'label': 'Variant', 'formatter': 'formatVariantString'}]
        colmodels = super().get_colmodels(remove_server_side_only=remove_server_side_only)
        return before_colmodels + colmodels


class UploadPipelineModifiedVariantsGrid(JqGridUserRowConfig):
    model = ModifiedImportedVariant
    caption = 'Modified Imported Variant'
    fields = ["variant__variantannotation__transcript_version__gene_version__gene_id", "variant__variantannotation__transcript_version__gene_version__gene_symbol__symbol",
              'old_multiallelic', 'old_variant']

    colmodel_overrides = {"variant__variantannotation__transcript_version__gene_version__gene_id": {"hidden": True},
                          "variant__variantannotation__transcript_version__gene_version__gene_symbol__symbol": {'label': 'Gene', 'formatter': 'geneLinkFormatter'}}

    def __init__(self, user, upload_pipeline_id):
        super().__init__(user)
        upload_pipeline = get_object_or_404(UploadPipeline, pk=upload_pipeline_id)

        queryset = get_queryset_for_latest_annotation_version(self.model, upload_pipeline.genome_build)
        queryset = queryset.filter(import_info__upload_step__upload_pipeline=upload_pipeline)
        queryset = Variant.annotate_variant_string(queryset, path_to_variant="variant__")

        field_names = self.get_field_names() + ["variant_string"]
        self.queryset = queryset.values(*field_names)
        self.extra_config.update({'sortname': 'variant_string',
                                  'sortorder': 'asc'})

    def get_colmodels(self, remove_server_side_only=False):
        before_colmodels = [{'index': 'variant_string', 'name': 'variant_string', 'label': 'Variant'}]
        colmodels = super().get_colmodels(remove_server_side_only=remove_server_side_only)
        return before_colmodels + colmodels
