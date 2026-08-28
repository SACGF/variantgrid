from functools import partial
from typing import Any

from django.db.models import QuerySet
from django.http import HttpRequest
from django.shortcuts import get_object_or_404

from annotation.annotation_version_querysets import get_queryset_for_latest_annotation_version
from snpdb.grids import AbstractSkippedAnnotationColumns
from snpdb.models import GenomeBuild, ProcessingStatus
from snpdb.models.models_variant import Variant
from snpdb.views.datatable_view import CellData, DatatableConfig, RichColumn, SortOrder
from upload.models import ModifiedImportedVariant, UploadPipeline, UploadStep, VCFPipelineStage


class UploadStepColumns(DatatableConfig[UploadStep]):

    def get_initial_queryset(self) -> QuerySet[UploadStep]:
        upload_pipeline_id = self.get_query_param("upload_pipeline")
        upload_pipeline = get_object_or_404(UploadPipeline, pk=upload_pipeline_id)
        upload_pipeline.file_upload.check_can_view(self.user)
        return UploadStep.objects.filter(upload_pipeline_id=upload_pipeline_id)

    @staticmethod
    def render_status(row: CellData):
        return ProcessingStatus(row["status"]).label

    @staticmethod
    def _render_pipeline_stage(column_name, row: CellData):
        value = None
        if cell := row[column_name]:
            value = VCFPipelineStage(cell).label
        return value

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
            RichColumn(key='pipeline_stage', orderable=True, renderer=partial(self._render_pipeline_stage, 'pipeline_stage')),
            RichColumn(key='pipeline_stage_dependency', orderable=True, renderer=partial(self._render_pipeline_stage, 'pipeline_stage_dependency')),
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


class UploadPipelineSkippedAnnotationColumns(AbstractSkippedAnnotationColumns):
    def _get_variant_source(self) -> tuple[Any, GenomeBuild]:
        upload_pipeline = self._get_upload_pipeline()
        return upload_pipeline.uploadedvcf.vcf, upload_pipeline.genome_build

    def _get_upload_pipeline(self) -> UploadPipeline:
        return get_object_or_404(UploadPipeline, pk=self.get_query_param("upload_pipeline_id"))


class UploadPipelineModifiedVariantsColumns(DatatableConfig[ModifiedImportedVariant]):
    """ Variants changed by decompose/normalise during import """
    grid_name = "Modified Imported Variant"
    GENE_SYMBOL_PATH = "variant__variantannotation__transcript_version__gene_version__gene_symbol__symbol"

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="variant_string", label="Variant", orderable=True, default_sort=SortOrder.ASC),
            RichColumn(key="operation", label="Operation", orderable=True),
            RichColumn(key=self.GENE_SYMBOL_PATH, name="gene_symbol", label="Gene", orderable=True,
                       client_renderer='renderGeneSymbol'),
            RichColumn(key="old_multiallelic", label="Old Multiallelic", orderable=True),
            RichColumn(key="old_variant", label="Old Variant", orderable=True),
            RichColumn(key="operation_detail", label="Operation Detail", orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[ModifiedImportedVariant]:
        upload_pipeline = get_object_or_404(UploadPipeline, pk=self.get_query_param("upload_pipeline_id"))
        qs = get_queryset_for_latest_annotation_version(ModifiedImportedVariant, upload_pipeline.genome_build)
        qs = qs.filter(import_info__upload_step__upload_pipeline=upload_pipeline)
        return Variant.annotate_variant_string(qs, path_to_variant="variant__")
