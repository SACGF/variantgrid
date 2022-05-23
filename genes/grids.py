import operator
from functools import reduce
from typing import Dict, Any

from django.conf import settings
from django.contrib.postgres.aggregates.general import StringAgg
from django.core.exceptions import PermissionDenied
from django.db.models import Count, TextField, QuerySet
from django.http import HttpRequest
from django.shortcuts import get_object_or_404
from django.urls.base import reverse

from analysis.models import VariantTag
from annotation.models.models import AnnotationVersion, GeneAnnotationVersion, InvalidAnnotationVersionError
from genes.models import CanonicalTranscript, GeneListCategory, GeneList, GeneSymbol, \
    GeneCoverageCanonicalTranscript, CanonicalTranscriptCollection, GeneCoverageCollection, TranscriptVersion, \
    GeneListGeneSymbol, GeneAnnotationRelease, ReleaseGeneVersion, GeneSymbolWiki
from library.django_utils.jqgrid_view import JQGridViewOp
from library.jqgrid_user_row_config import JqGridUserRowConfig
from snpdb.grid_columns.custom_columns import get_custom_column_fields_override_and_sample_position
from snpdb.grids import AbstractVariantGrid
from snpdb.models import UserSettings, Q, VariantGridColumn, Tag
from snpdb.models.models_genome import GenomeBuild
from snpdb.variant_queries import get_variant_queryset_for_gene_symbol, variant_qs_filter_has_internal_data
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder
from uicore.json.json_types import JsonDataType


class GeneListsGrid(JqGridUserRowConfig):
    model = GeneList
    caption = 'Gene Lists'
    # Category is only shown when gene_id provided (hidden set in get_colmodels)
    fields = ["id", "category__name", "name", "user__username", "import_status"]
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {"width": 400,
                 'formatter': 'linkFormatter',
                 'formatter_kwargs': {"icon_css_class": "gene-list-icon",
                                      "url_name": "view_gene_list",
                                      "url_object_column": "id"}},
        'user__username': {'label': 'Uploaded by'},
        'category__name': {"label": "Category", "hidden": True}
    }

    def __init__(self, user, **kwargs):
        super().__init__(user)

        queryset = GeneList.filter_for_user(user, success_only=False)
        gene_list_category_id = kwargs.pop("gene_list_category_id", None)
        gene_symbol_id = kwargs.pop("gene_symbol", None)
        self.show_category = bool(gene_symbol_id)

        if gene_symbol_id:
            gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol_id)
            queryset = GeneList.visible_gene_lists_containing_gene_symbol(queryset, gene_symbol)
        else:
            if gene_list_category_id:
                gene_list_category_id = int(gene_list_category_id)
                gene_list_category = GeneListCategory.objects.get(pk=gene_list_category_id)
                if gene_list_category.public:
                    queryset = GeneList.objects.all()
                queryset = queryset.filter(category=gene_list_category)
            else:
                queryset = queryset.filter(category__isnull=True)

        queryset = queryset.annotate(num_genes=Count("genelistgenesymbol"))
        field_names = self.get_field_names() + ["num_genes"]
        self.queryset = queryset.values(*field_names)

    def get_colmodels(self, remove_server_side_only=False):
        colmodels = super().get_colmodels(remove_server_side_only=remove_server_side_only)
        gene_list_genes_colmodel = {'index': 'num_genes', 'name': 'num_genes', 'label': '# gene symbols'}
        colmodels += [gene_list_genes_colmodel]
        if self.show_category:
            for cm in colmodels:
                if cm['index'] == 'category__name':
                    cm['hidden'] = False
        return colmodels


class GeneListGenesGrid(JqGridUserRowConfig):
    model = GeneListGeneSymbol
    caption = 'Gene List Genes'
    colmodel_overrides = {
        'gene_symbol_alias__source': {"label": "Alias"},
    }
    fields = ["original_name", "gene_symbol__symbol", "gene_symbol_alias__source"]

    def __init__(self, user, gene_list_id):
        super().__init__(user)
        gene_list = GeneList.get_for_user(user, gene_list_id, success_only=False)
        queryset = self.model.objects.filter(gene_list=gene_list)

        annotation_kwargs = {}
        self.annotation_field_labels = {}
        for release in GeneAnnotationRelease.get_for_latest_annotation_versions_for_builds():
            field_name = f"release_{release.pk}"
            self.annotation_field_labels[field_name] = str(release)
            annotation_kwargs[field_name] = GeneListGeneSymbol.get_joined_genes_qs_annotation_for_release(release)
        queryset = queryset.annotate(**annotation_kwargs).order_by("original_name")
        field_names = self.get_field_names() + list(sorted(self.annotation_field_labels))
        self.queryset = queryset.values(*field_names)
        self.extra_config.update({'sortname': 'original_name',
                                  'sortorder': 'asc'})

    def get_colmodels(self, remove_server_side_only=False):
        colmodels = super().get_colmodels(remove_server_side_only=False)
        for field_name, label in self.annotation_field_labels.items():
            cm = {'index': field_name, 'name': field_name, 'label': label, "width": 400}
            colmodels.append(cm)
        return colmodels


class GeneSymbolVariantsGrid(AbstractVariantGrid):
    """ Uses custom columns subtracting away the gene annotations (as they're displayed above) """
    caption = 'Gene Variants'
    fields = ["id", "locus__contig__name", 'locus__position', 'locus__ref', 'alt']
    colmodel_overrides = {
        'id': {'editable': False, 'width': 90, 'fixed': True, 'formatter': 'detailsLink'},
        'tags_global': {'classes': 'no-word-wrap', 'formatter': 'tagsGlobalFormatter', 'sortable': False},
    }

    def __init__(self, user, gene_symbol, genome_build_name, **kwargs):
        extra_filters = kwargs.pop("extra_filters", None)
        super().__init__(user)

        user_settings = UserSettings.get_for_user(user)
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        annotation_version = AnnotationVersion.latest(genome_build)
        fields, override, _ = get_custom_column_fields_override_and_sample_position(user_settings.columns,
                                                                                    annotation_version)
        self.fields = self._get_non_gene_fields(fields)
        self.update_overrides(override)

        gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol)
        genes_qs = get_variant_queryset_for_gene_symbol(gene_symbol, annotation_version)
        queryset = variant_qs_filter_has_internal_data(genes_qs, genome_build)
        if extra_filters:
            # Hotspot filters
            if protein_position := extra_filters.get("protein_position"):
                transcript_version = TranscriptVersion.objects.get(pk=extra_filters["protein_position_transcript_version_id"])
                queryset = queryset.filter(varianttranscriptannotation__transcript_version=transcript_version,
                                           varianttranscriptannotation__protein_position__icontains=protein_position)
            tag_id = extra_filters.get("tag")
            if tag_id is not None:  # "" for all tags
                tags_qs = VariantTag.objects.all()
                if tag_id:
                    tag = get_object_or_404(Tag, pk=tag_id)
                    tags_qs = tags_qs.filter(tag=tag)
                alleles_qs = tags_qs.values_list("variant__variantallele__allele")
                queryset = queryset.filter(variantallele__allele__in=alleles_qs)

        self.queryset = queryset.distinct().values(*self.get_queryset_field_names())
        self.extra_config.update({'sortname': "locus__position",
                                  'sortorder': "asc",
                                  'shrinkToFit': False})

    @staticmethod
    def _get_non_gene_fields(fields):
        """ Remove fields that'll all be the same """
        non_gene_fields = []
        for f in fields:
            if f == "tags":
                continue
            keep = True
            for gene_fields in ["__transcript_version__", "__gene__"]:
                if gene_fields in f:
                    keep = False
                    break
            if keep:
                non_gene_fields.append(f)
        return non_gene_fields


def _get_gene_fields():
    q_gene = Q(variant_column__contains='__gene__') | Q(variant_column__contains='__gene_version__')
    columns_qs = VariantGridColumn.objects.filter(q_gene).order_by("pk")
    first_fields = ["gene_version__gene_symbol__symbol", "gene_version__gene__identifier", "gene_version__version"]
    fields = []
    for variant_column in columns_qs.values_list("variant_column", flat=True):
        gene_column = variant_column.replace("variantannotation__", "").replace("transcript_version__", "")
        if gene_column.startswith("gene__"):
            gene_column = "gene_version__" + gene_column
        if gene_column not in first_fields:
            fields.append(gene_column)

    return first_fields + fields


class GenesGrid(JqGridUserRowConfig):
    model = ReleaseGeneVersion
    caption = "Gene Release"
    colmodel_overrides = {
        'gene_version__gene_symbol__symbol': {'formatter': 'geneSymbolLink'},
        "gene_version__hgnc__gene_symbol__symbol": {"label": "HGNC Symbol"},
    }

    def __init__(self, user, genome_build_name, **kwargs):
        extra_filters = kwargs.pop("extra_filters", None)
        self.fields = _get_gene_fields()
        super().__init__(user)
        queryset = self.model.objects.all()
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

        # TODO: Put back all the other fields - joining to HGNC and to GeneAnnotation

        av = AnnotationVersion.latest(genome_build)
        gene_annotation_release_id = av.gene_annotation_version.gene_annotation_release_id
        if extra_filters:
            gene_annotation_release_id = extra_filters["gene_annotation_release_id"]
            if column := extra_filters.get("column"):
                if column in self.fields:
                    is_null = extra_filters.get("is_null", False)
                    kwargs = {f"{column}__isnull": is_null}
                    queryset = queryset.filter(**kwargs)
                else:
                    raise PermissionDenied(f"Bad column '{column}'")

        gene_annotation_version = GeneAnnotationVersion.objects.filter(gene_annotation_release_id=gene_annotation_release_id).order_by("annotation_date").last()
        if gene_annotation_version is None:
            raise InvalidAnnotationVersionError(f"No gene annotation version for gene_annotation_release: {gene_annotation_release_id}")
        queryset = queryset.filter(release_id=gene_annotation_release_id)
        queryset = queryset.filter(gene_version__gene__geneannotation__version=gene_annotation_version)
        self.queryset = queryset.values(*self.get_field_names())
        grid_export_url = reverse("genes_grid", kwargs={"genome_build_name": genome_build_name,
                                                        "op": JQGridViewOp.DOWNLOAD})
        self.extra_config.update({'sortname': 'gene_version__gene_symbol',
                                  'sortorder': 'asc',
                                  'shrinkToFit': False,
                                  'grid_export_url': grid_export_url})


class CanonicalTranscriptCollectionsGrid(JqGridUserRowConfig):
    model = CanonicalTranscriptCollection
    caption = 'CanonicalTranscripts'
    fields = ["id", "description", "annotation_consortium"]
    colmodel_overrides = {'id': {'editable': False, 'width': 90, 'fixed': True,
                                 'formatter': 'viewCanonicalTranscriptCollection'}}

    def __init__(self, user):
        super().__init__(user)

        queryset = self.model.objects.all()
        queryset = queryset.annotate(enrichment_kits=StringAgg("enrichmentkit__name", ',', output_field=TextField()))
        field_names = self.get_field_names() + ["enrichment_kits"]
        self.queryset = queryset.values(*field_names)

    def get_colmodels(self, remove_server_side_only=False):
        colmodels = super().get_colmodels(remove_server_side_only=remove_server_side_only)
        enrichment_kits_colmodel = {'index': 'enrichment_kits', 'name': 'enrichment_kits',
                                    'label': 'Enrichment Kits', 'width': 230}
        colmodels += [enrichment_kits_colmodel]
        return colmodels


class CanonicalTranscriptGrid(JqGridUserRowConfig):
    model = CanonicalTranscript
    caption = 'CanonicalTranscripts'
    fields = ["gene_symbol__symbol", "transcript__identifier", "original_gene_symbol", "original_transcript"]
    colmodel_overrides = {'gene_symbol__symbol': {'label': "Matched Symbol"},
                          'transcript__identifier': {'label': "Matched Transcript"}}

    def __init__(self, user, pk):
        super().__init__(user)
        canonical_transcript_collection = get_object_or_404(CanonicalTranscriptCollection, pk=pk)
        queryset = self.model.objects.all()
        queryset = queryset.filter(collection=canonical_transcript_collection)

        self.queryset = queryset.values(*self.get_field_names())


class QCGeneCoverageGrid(JqGridUserRowConfig):
    model = GeneCoverageCanonicalTranscript
    caption = 'QC'
    fields = ["gene_symbol__symbol", "transcript__identifier", "original_gene_symbol", "original_transcript", "min",
              "mean", "std_dev", "percent_1x", "percent_10x", "percent_20x"]
    number_format = {'formatter': 'number', 'width': 80}
    colmodel_overrides = {'gene_symbol__symbol': {"width": 110},
                          'transcript__identifier': {"width": 110},
                          "original_gene_symbol": {'label': 'original symbol'},
                          "original_transcript": {'label': 'original transcript'},
                          'min': {'width': 40},
                          'mean': number_format,
                          'std_dev': number_format,
                          'percent_1x': number_format,
                          'percent_10x': number_format,
                          'percent_20x': number_format}

    def __init__(self, user, gene_coverage_collection_id, gene_list_id_list=None):
        super().__init__(user)
        gene_coverage_collection = get_object_or_404(GeneCoverageCollection, pk=gene_coverage_collection_id)
        gene_symbols = set()
        if gene_list_id_list:
            gene_list_ids = gene_list_id_list.split("/")
            if gene_list_ids:
                for gene_list_id in gene_list_ids:
                    gene_list = get_object_or_404(GeneList, pk=gene_list_id)
                    gene_symbols.update(gene_list.get_gene_names())

        q = self.get_coverage_q(gene_coverage_collection, gene_symbols)
        queryset = self.model.objects.filter(q)
        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': 'id',
                                  'sortorder': 'desc'})

    def get_coverage_q(self, gene_coverage_collection, gene_symbols) -> Q:
        filters = [
            Q(gene_coverage_collection=gene_coverage_collection),
        ]
        if gene_symbols:
            filters.append(Q(gene_symbol__in=gene_symbols))

        return reduce(operator.and_, filters)


class UncoveredGenesGrid(QCGeneCoverageGrid):

    def __init__(self, user, min_depth=settings.SEQAUTO_MIN_COVERAGE, **kwargs):
        self.min_depth = min_depth
        super().__init__(user, **kwargs)

    def get_coverage_q(self, gene_coverage_collection, gene_symbols) -> Q:
        q = super().get_coverage_q(gene_coverage_collection, gene_symbols)
        return q & Q(min__lt=self.min_depth)


class GeneSymbolWikiColumns(DatatableConfig[GeneSymbolWiki]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.download_csv_button_enabled = True

        self.rich_columns = [
            RichColumn('gene_symbol', renderer=self.render_gene_symbol, client_renderer="renderGeneSymbol"),
            RichColumn('markdown'),
            RichColumn('last_edited_by__username', name='user', orderable=True),
            RichColumn('created', client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn('modified', client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
        ]

    @staticmethod
    def render_gene_symbol(row: Dict[str, Any]) -> JsonDataType:
        gene_symbol = row["gene_symbol"]
        return {"id": gene_symbol}

    def get_initial_queryset(self) -> QuerySet[GeneSymbolWiki]:
        return GeneSymbolWiki.objects.all()
