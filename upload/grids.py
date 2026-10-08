from functools import partial
from typing import Any, Optional

from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.http import HttpRequest
from django.shortcuts import get_object_or_404
from django.urls import reverse

from annotation.annotation_version_querysets import get_queryset_for_latest_annotation_version
from snpdb.grids import AbstractSkippedAnnotationColumns
from snpdb.models import GenomeBuild, ProcessingStatus
from snpdb.models.models_variant import Variant
from snpdb.views.datatable_view import CellData, DatatableConfig, RichColumn, SortOrder
from upload.file_type_icons import file_type_icon_html
from upload.models import (
    FileUpload,
    ModifiedImportedVariant,
    UploadedFileTypes,
    UploadedVCF,
    UploadPipeline,
    UploadStep,
    VCFPipelineStage,
)
from upload.uploaded_file_type import get_upload_data_by_file_upload_id
from upload.views.views_json import get_remaining_annotation_runs


def get_upload_pipeline_for_user(user: User, upload_pipeline_id) -> UploadPipeline:
    upload_pipeline = get_object_or_404(UploadPipeline, pk=upload_pipeline_id)
    upload_pipeline.file_upload.check_can_view(user)
    return upload_pipeline


def get_status_icon(status, requires_user_input_url: Optional[str] = None) -> dict:
    """ requires_user_input_url: the data's page, where the user sets what the import is waiting on """
    if requires_user_input_url:
        return {'icon': 'fa-exclamation-triangle', 'css': 'text-warning',
                'title': 'Requires input - set genome build', 'url': requires_user_input_url}

    ICONS = {
        ProcessingStatus.CREATED: {'icon': 'fa-clock', 'title': 'Queued'},
        ProcessingStatus.PROCESSING: {'icon': 'fa-spinner fa-spin', 'title': 'Processing'},
        ProcessingStatus.ERROR: {'icon': 'fa-times-circle', 'css': 'text-danger', 'title': 'Error'},
        ProcessingStatus.SUCCESS: {'icon': 'fa-check-circle', 'css': 'text-success', 'title': 'Success'},
        ProcessingStatus.TERMINATED_EARLY: {'icon': 'fa-exclamation-triangle', 'css': 'text-warning', 'title': 'Terminated early'},
    }
    return ICONS.get(status, {})


class FileUploadColumns(DatatableConfig[FileUpload]):
    """ The upload page's files - the user's own, or everyone's for a superuser """
    grid_name = "File Uploads"
    search_box_enabled = True
    # Files the system makes for itself - shown only when picked from the page's file type select
    INTERNAL_FILE_TYPES = {UploadedFileTypes.LIFTOVER, UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY,
                           UploadedFileTypes.GENE_LEVEL_INSERT_VARIANTS_ONLY}
    ALL_FILE_TYPES = "all"

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self._upload_data_by_file_upload_id = {}
        self.rich_columns = [
            RichColumn(key="id", visible=False, search=False),
            RichColumn(key="uploadpipeline__status", name="status", label="Status", orderable=True, search=False,
                       extra_columns=["id", "file_type"], renderer=self.render_status,
                       client_renderer="renderUploadStatus"),
            RichColumn(key="name", label="Name", orderable=True, extra_columns=["id", "uploadpipeline__id"],
                       renderer=self.render_name, client_renderer="TableFormat.linkUrl"),
            RichColumn(key="file_type", label="File Type", orderable=True, search=False,
                       extra_columns=["id", "uploadpipeline__status"], renderer=self.render_file_type,
                       client_renderer="renderUploadFileType"),
            RichColumn(name="genome_build", label="Genome Build", extra_columns=["id", "uploadpipeline__status"],
                       renderer=self.render_genome_build),
            self.user_column(label="User", enabled=self.user.is_superuser),
            RichColumn(key="created", label="Uploaded", orderable=True, search=False,
                       default_sort=SortOrder.DESC, client_renderer="TableFormat.timeAgo"),
            RichColumn(key="id", name="delete", label="", search=False, include_in_csv=False,
                       renderer=self.render_delete_url, client_renderer="TableFormat.deleteRow"),
        ]

    def get_initial_queryset(self) -> QuerySet[FileUpload]:
        qs = FileUpload.objects.all()
        if not self.user.is_superuser:
            qs = qs.filter(user=self.user)
        return qs

    def filter_queryset(self, qs: QuerySet[FileUpload]) -> QuerySet[FileUpload]:
        file_type = self.get_query_param("file_type")
        if not file_type:
            qs = qs.exclude(file_type__in=self.INTERNAL_FILE_TYPES)
        elif file_type != self.ALL_FILE_TYPES:
            qs = qs.filter(file_type=file_type)
        return qs

    def pre_render(self, qs: QuerySet[FileUpload], rows: list[dict]):
        super().pre_render(qs, rows)
        file_type_by_file_upload_id = {row["id"]: row["file_type"] for row in rows}
        self._upload_data_by_file_upload_id = get_upload_data_by_file_upload_id(file_type_by_file_upload_id)

    def render_status(self, row: CellData) -> dict:
        status = row["uploadpipeline__status"] or ProcessingStatus.ERROR
        requires_user_input_url = None
        upload_data = self._upload_data_by_file_upload_id.get(row["id"])
        if upload_data and upload_data.requires_user_input:
            requires_user_input_url = upload_data.get_data_url()
        status_icon = get_status_icon(status, requires_user_input_url)
        if not row["file_type"]:
            status_icon["title"] = "Could not determine how to read file"
        return {"status": status, **status_icon}

    @staticmethod
    def render_name(row: CellData) -> dict:
        if upload_pipeline_id := row["uploadpipeline__id"]:
            url = reverse('view_upload_pipeline', kwargs={'upload_pipeline_id': upload_pipeline_id})
        else:
            url = reverse('view_uploaded_file', kwargs={'file_upload_id': row["id"]})
        return {"text": row["name"], "url": url}

    def render_file_type(self, row: CellData) -> dict:
        file_type = row["file_type"]
        data = {"icon": file_type_icon_html(file_type)}
        if file_type:
            data["label"] = UploadedFileTypes(file_type).label
        if row["uploadpipeline__status"] in (ProcessingStatus.SUCCESS, ProcessingStatus.TERMINATED_EARLY):
            if upload_data := self._upload_data_by_file_upload_id.get(row["id"]):
                data["url"] = upload_data.get_data_url()
        return data

    def render_genome_build(self, row: CellData) -> str:
        upload_data = self._upload_data_by_file_upload_id.get(row["id"])
        if not (upload_data and (genome_build := upload_data.genome_build)):
            return ""
        text = str(genome_build)
        if row["uploadpipeline__status"] == ProcessingStatus.PROCESSING and isinstance(upload_data, UploadedVCF):
            if remaining_annotation_runs := get_remaining_annotation_runs(upload_data, genome_build):
                text += f" (annotating: {remaining_annotation_runs} runs remaining)"
        return text

    @staticmethod
    def render_delete_url(row: CellData) -> str:
        return reverse('upload_file_delete', kwargs={'pk': row.value})


class UploadStepColumns(DatatableConfig[UploadStep]):

    def get_initial_queryset(self) -> QuerySet[UploadStep]:
        upload_pipeline = get_upload_pipeline_for_user(self.user, self.get_query_param("upload_pipeline"))
        return UploadStep.objects.filter(upload_pipeline=upload_pipeline)

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
        return get_upload_pipeline_for_user(self.user, self.get_query_param("upload_pipeline_id"))


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
        upload_pipeline = get_upload_pipeline_for_user(self.user, self.get_query_param("upload_pipeline_id"))
        qs = get_queryset_for_latest_annotation_version(ModifiedImportedVariant, upload_pipeline.genome_build)
        qs = qs.filter(import_info__upload_step__upload_pipeline=upload_pipeline)
        return Variant.annotate_variant_string(qs, path_to_variant="variant__")
