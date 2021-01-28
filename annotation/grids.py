from django.shortcuts import get_object_or_404

from annotation.models import VariantAnnotationVersion, AnnotationRun, HumanProteinAtlasAbundance
from annotation.models.models import HumanProteinAtlasAnnotationVersion, HumanProteinAtlasTissueSample
from genes.models import GeneVersion
from library.jqgrid_abstract_genes_grid import AbstractGenesGrid
from library.jqgrid_user_row_config import JqGridUserRowConfig
from snpdb.models.models_genome import GenomeBuild


class AnnotationRunGrid(JqGridUserRowConfig):
    model = AnnotationRun
    caption = 'AnnotationRun'
    fields = [
        'id', 'status',
        'annotation_range_lock__version__genome_build__name',
        'annotation_range_lock__version', 'annotation_range_lock__count', 'vep_skipped_count',
        'annotation_range_lock__min_variant', 'annotation_range_lock__max_variant',
        'task_id', 'created', 'dump_start', 'dump_end',
        'annotation_start', 'annotation_end', 'upload_start', 'upload_end',
        'error_exception', 'vcf_dump_filename', 'vcf_annotated_filename'
    ]  # 'pipeline_stdout' and 'pipeline_stderr' are too long!
    colmodel_overrides = {
        'id': {'width': 20, 'formatter': 'viewAnnotationRunLink'},
        'annotation_range_lock__version__genome_build__name': {'label': 'Build'},
        'annotation_range_lock__version': {'label': 'Version'},
        'annotation_range_lock__count': {'label': 'Variant Count'},
        'annotation_range_lock__min_variant': {'label': 'Min Variant'},
        'annotation_range_lock__max_variant': {'label': 'Max Variant'}
    }

    def __init__(self, user, genome_build_name):
        super().__init__(user)
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        queryset = self.model.objects.filter(annotation_range_lock__version__genome_build=genome_build)
        self.queryset = queryset.values(*self.get_field_names())

        self.extra_config.update({'sortname': "id",
                                  'sortorder': "desc",
                                  'shrinkToFit': False})


class VariantAnnotationVersionGrid(JqGridUserRowConfig):
    model = VariantAnnotationVersion
    caption = 'VariantAnnotationVersion'
    fields = [
        'id', "vep", "annotation_consortium", 'created', 'last_checked_date', 'gene_annotation_release',
        "ensembl", "ensembl_funcgen", "ensembl_variation", "ensembl_io",
        "thousand_genomes", "cosmic", "hgmd", "assembly", "dbsnp",
        "gencode", "genebuild", "gnomad", "refseq", "regbuild", "sift", "dbnsfp"
    ]

    def __init__(self, user, genome_build_name):
        super().__init__(user)
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        queryset = self.model.objects.filter(genome_build=genome_build)
        self.queryset = queryset.values(*self.get_field_names())

        self.extra_config.update({'sortname': "created",
                                  'sortorder': "desc",
                                  'shrinkToFit': False})


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
