from django.db.models import QuerySet
from django.http import HttpRequest
from django.shortcuts import get_object_or_404

from annotation.annotation_version_querysets import get_queryset_for_latest_annotation_version
from library.jqgrid.jqgrid_user_row_config import JqGridUserRowConfig
from snpdb.models import ProcessingStatus
from snpdb.models.models_variant import Variant
from snpdb.views.datatable_view import DatatableConfig, RichColumn, CellData
from upload.models import UploadStep, ModifiedImportedVariant, UploadPipeline


class UploadStepColumns(DatatableConfig[UploadStep]):

    def get_initial_queryset(self) -> QuerySet[UploadStep]:
        upload_pipeline_id = self.get_query_param("upload_pipeline")
        upload_pipeline = get_object_or_404(UploadPipeline, pk=upload_pipeline_id)
        upload_pipeline.uploaded_file.check_can_view(self.user)
        return UploadStep.objects.filter(upload_pipeline_id=upload_pipeline_id)

    @staticmethod
    def render_status(row: CellData):
        return ProcessingStatus(row["status"]).label

    @staticmethod
    def render_duration(row: CellData):
        start_date = row["start_date"]
        end_date = row["end_date"]
        if start_date and end_date:
            delta = end_date - start_date
            return f"{delta.total_seconds():.2f}"
        else:
            return ""

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.scroll_x = True
        self.expand_client_renderer = DatatableConfig._row_expand_ajax('upload_step_detail',
                                                                       expected_height=120)
        self.rich_columns = [
            RichColumn(key='sort_order', orderable=True),
            RichColumn(key='id', orderable=True),
            RichColumn(key='name', orderable=True),
            RichColumn(key='status', orderable=True, renderer=UploadStepColumns.render_status),
            RichColumn(key='items_processed', css_class='num', orderable=True),
            RichColumn(key='error_message', orderable=True),
            RichColumn(key='input_filename', orderable=True),
            RichColumn(key='output_filename', orderable=True),
            RichColumn(key='start_date', client_renderer='TableFormat.timestampMilliseconds', orderable=True),
            RichColumn(key='end_date', client_renderer='TableFormat.timestampMilliseconds', orderable=True),
            RichColumn(name='duration', label="Duration Seconds", extra_columns=["start_date", "end_date"],
                       renderer=UploadStepColumns.render_duration, css_class="num")
        ]


class UploadPipelineSkippedAnnotationGrid(JqGridUserRowConfig):
    model = Variant
    caption = 'Skipped Annotation'
    fields = ["id", "variantannotation__vep_skipped_reason", "variantannotation__annotation_run_id"]

    colmodel_overrides = {"id": {"hidden": True},
                          "variantannotation__annotation_run_id": {'formatter': 'formatAnnotationRunLink'}}

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
    fields = ["variant__variantannotation__transcript_version__gene_version__gene__identifier", "variant__variantannotation__transcript_version__gene_version__gene_symbol__symbol",
              'old_multiallelic', 'old_variant']

    colmodel_overrides = {"variant__variantannotation__transcript_version__gene_version__gene__identifier": {"hidden": True},
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
