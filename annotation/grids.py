import json
from datetime import timedelta
from typing import Dict, Any

from django.db.models import QuerySet, ExpressionWrapper, F, fields

from annotation.models import VariantAnnotationVersion, AnnotationRun, AnnotationStatus
from genes.models_enums import AnnotationConsortium
from library.jqgrid_user_row_config import JqGridUserRowConfig
from snpdb.models.models_genome import GenomeBuild
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder, CellData


class AnnotationRunColumns(DatatableConfig):

    def status(self, row: Dict[str, Any]):
        return AnnotationStatus(row["status"]).label

    def format_timedelta(self, cell: CellData):
        delta: timedelta = cell.value
        if delta is None:
            return '-'
        seconds = delta.total_seconds()
        hours, remainder = divmod(seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        if int(seconds) == 0:
            return '< 1 second'
        return '{:02}:{:02}:{:02}'.format(int(hours), int(minutes), int(seconds))

    def __init__(self, request):
        super().__init__(request)

        preview_columns = ["error_exception", "pipeline_stdout", "pipeline_stderr"]

        self.rich_columns = [
            RichColumn(key="id", label='ID', orderable=True, client_renderer='idRenderer'),
            RichColumn(key="status", orderable=True, renderer=self.status),
            RichColumn(key="annotation_range_lock__version__genome_build__name", label='Build', orderable=True),
            RichColumn(key="annotation_range_lock__version__id", label='Version', orderable=True),
            RichColumn(key="annotation_range_lock__count", label='Var Count', orderable=True),
            RichColumn(key="vep_skipped_count", label="VEP Skipped", orderable=True),
            RichColumn(key="annotation_range_lock__min_variant__id", label="Min Var", orderable=True),
            RichColumn(key="annotation_range_lock__max_variant__id", label="Max Var", orderable=True),

            # RichColumn(key="task_id", label="Task ID", orderable=True),
            RichColumn(key="created", client_renderer='TableFormat.timestamp', orderable=True),
            # RichColumn(key="dump_start", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="dump_duration", renderer=self.format_timedelta, orderable=True),
            # RichColumn(key="annotation_start", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="annotation_duration", renderer=self.format_timedelta, orderable=True),
            # RichColumn(key="upload_start", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="upload_duration", renderer=self.format_timedelta, orderable=True),

            RichColumn(key=None, name='preview', label='Data',
                       renderer=lambda x: {key: x[key] for key in preview_columns},
                       client_renderer=f'TableFormat.preview.bind(null, {json.dumps(preview_columns)})',
                       extra_columns=preview_columns)
        ]

    def get_initial_queryset(self):
        qs = AnnotationRun.objects.all()
        qs = qs.annotate(
            dump_duration=ExpressionWrapper(F('dump_end') - F('dump_start'), output_field=fields.DurationField()))
        qs = qs.annotate(
            annotation_duration=ExpressionWrapper(F('annotation_end') - F('annotation_start'), output_field=fields.DurationField()))
        qs = qs.annotate(
            upload_duration=ExpressionWrapper(F('upload_end') - F('upload_start'), output_field=fields.DurationField()))
        return qs

    def filter_queryset(self, qs: QuerySet) -> QuerySet:
        if variant_annotation_version_id := self.get_query_param("variant_annotation_version_id"):
            variant_annotation_version = VariantAnnotationVersion.objects.get(pk=variant_annotation_version_id)
            qs = qs.filter(annotation_range_lock__version=variant_annotation_version)
        if status_str := self.get_query_param("status"):
            if status_str == "outstanding":
                qs = qs.exclude(status__in={AnnotationStatus.FINISHED})
        return qs


class VariantAnnotationVersionGrid(JqGridUserRowConfig):
    model = VariantAnnotationVersion
    caption = 'VariantAnnotationVersion'
    fields = [
        'id', "vep", "annotation_consortium", 'created', 'last_checked_date', 'gene_annotation_release__version',
        "ensembl", "ensembl_funcgen", "ensembl_variation", "ensembl_io",
        "thousand_genomes", "cosmic", "hgmd", "assembly", "dbsnp",
        "gencode", "genebuild", "gnomad", "refseq", "regbuild", "sift", "dbnsfp", "distance"
    ]
    colmodel_overrides = {
        'gene_annotation_release__version': {"label": "Gene Annotation Release"},
    }

    def __init__(self, user, genome_build_name):
        super().__init__(user)
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        queryset = self.model.objects.filter(genome_build=genome_build)
        self.queryset = queryset.values(*self.get_field_names())

        self.extra_config.update({'sortname': "created",
                                  'sortorder': "desc",
                                  'shrinkToFit': False})


class VariantAnnotationVersionColumns(DatatableConfig[VariantAnnotationVersion]):

    def __init__(self, request):
        super().__init__(request)

        self.rich_columns = [
            RichColumn(key="id", label="ID", orderable=True, default_sort=SortOrder.DESC, css_class='toggle-link'),
            RichColumn(key="vep", label="VEP", orderable=True),
            RichColumn(key="annotation_consortium", orderable=True, renderer=lambda x: AnnotationConsortium(x['annotation_consortium']).label),
            RichColumn(key="created", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="last_checked_date", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="columns_version", label="Columns Version", orderable=True),
            RichColumn(key="gene_annotation_release__version", label="Gene Annotation Release Version", orderable=True),
            RichColumn(key="ensembl", orderable=True),
            RichColumn(key="ensembl_funcgen", orderable=True),
            RichColumn(key="ensembl_variation", orderable=True, detail=True),
            RichColumn(key="ensembl_io", label="Ensembl-io", orderable=True, detail=True),
            RichColumn(key="thousand_genomes", orderable=True, detail=True),
            RichColumn(key="cosmic", orderable=True, detail=True),
            RichColumn(key="hgmd", label="HGMD", orderable=True, detail=True),
            RichColumn(key="assembly", orderable=True, detail=True),
            RichColumn(key="gencode", label="GENCODE", orderable=True, detail=True),
            RichColumn(key="genebuild", label="GeneBuild", orderable=True, detail=True),
            RichColumn(key="gnomad", label="gnomAD", orderable=True, detail=True),
            RichColumn(key="refseq", orderable=True, detail=True),
            RichColumn(key="regbuild", label="Regulatory Build", orderable=True, detail=True),
            RichColumn(key="sift", orderable=True, detail=True),
            RichColumn(key="dbnsfp", label="dbNSFP", orderable=True, detail=True),
            RichColumn(key="distance", orderable=True, detail=True)
        ]

    def get_initial_queryset(self) -> QuerySet[VariantAnnotationVersion]:
        genome_build = GenomeBuild.get_name_or_alias(self.get_query_param("genome_build_name"))
        return VariantAnnotationVersion.objects.filter(genome_build=genome_build)
