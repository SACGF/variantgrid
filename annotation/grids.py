import json
from datetime import timedelta
from typing import Dict, Any

from django.db.models import QuerySet, ExpressionWrapper, F, fields
from django.shortcuts import get_object_or_404

from annotation.models import VariantAnnotationVersion, AnnotationRun, HumanProteinAtlasAbundance, AnnotationStatus
from annotation.models.models import HumanProteinAtlasAnnotationVersion, HumanProteinAtlasTissueSample
from genes.models import GeneVersion
from genes.models_enums import AnnotationConsortium
from library.jqgrid_abstract_genes_grid import AbstractGenesGrid
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
        if genome_build_str := self.get_query_param("genome_build"):
            genome_build = GenomeBuild.get_name_or_alias(genome_build_str)
            qs = qs.filter(annotation_range_lock__version__genome_build=genome_build)
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


class VaraintAnnotationVersionColumns(DatatableConfig[VariantAnnotationVersion]):

    def __init__(self, request):
        super().__init__(request)

        self.rich_columns = [
            RichColumn(key="id", label="ID", orderable=True, default_sort=SortOrder.DESC, css_class='toggle-link'),
            RichColumn(key="vep", label="VEP", orderable=True),
            RichColumn(key="annotation_consortium", orderable=True, renderer=lambda x: AnnotationConsortium(x['annotation_consortium']).label),
            RichColumn(key="created", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="last_checked_date", client_renderer='TableFormat.timestamp', orderable=True),
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
            RichColumn(key="gnomad", label="gnomAD", orderable=True,  detail=True),
            RichColumn(key="refseq", orderable=True, detail=True),
            RichColumn(key="regbuild", label="Regulatory Build", orderable=True, detail=True),
            RichColumn(key="sift", orderable=True, detail=True),
            RichColumn(key="dbnsfp", label="dbNSFP", orderable=True, detail=True),
            RichColumn(key="distance", orderable=True, detail=True)
        ]

    def get_initial_queryset(self) -> QuerySet[VariantAnnotationVersion]:
        genome_build = GenomeBuild.get_name_or_alias(self.get_query_param("genome_build_name"))
        return VariantAnnotationVersion.objects.filter(genome_build=genome_build)


class TissueGeneGrid(AbstractGenesGrid):
    model = GeneVersion
    caption = "Tissue Genes"

    def __init__(self, user, human_protein_atlas_version_id, tissue_sample_id, min_abundance):
        super().__init__(user)

        self.human_protein_atlas_version = get_object_or_404(HumanProteinAtlasAnnotationVersion, pk=human_protein_atlas_version_id)
        self.tissue_sample = get_object_or_404(HumanProteinAtlasTissueSample, pk=tissue_sample_id)
        self.min_abundance = min_abundance
        self.extra_config.update({'sortname': 'name',
                                  'sortorder': 'asc'})

    def get_column_names(self):
        return ["gene_id"]

    def get_sql_params_and_columns(self, request):
        sql_template = """  select distinct genes_geneversion.gene_symbol_id as name, annotation_humanproteinatlasannotation.gene_id as gene_id
                            from annotation_humanproteinatlasannotation
                            join genes_gene on (genes_gene.identifier=annotation_humanproteinatlasannotation.gene_id)
                            where annotation_humanproteinatlasannotation.tissue_sample_id = %s
                            AND annotation_humanproteinatlasannotation.abundance in %s
        """

        from_table = "annotation_humanproteinatlasannotation"
        to_table = self.human_protein_atlas_version.get_partition_table()
        sql = sql_template.replace(from_table, to_table)

        abundances = HumanProteinAtlasAbundance.get_abundance_or_greater_levels(self.min_abundance)
        params = [self.tissue_sample.pk, tuple(abundances)]
        return (sql, params, None, False)

    def get_labels(self):
        return ['Gene ID', "Ensembl Gene ID"]

    @property
    def csv_name(self):
        return f"{self.tissue_sample.name}_{self.min_abundance}_tissue_genes"
