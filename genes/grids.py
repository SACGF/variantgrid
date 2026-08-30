import operator
import os
from functools import reduce
from typing import Any, Optional

from django.conf import settings
from django.core.exceptions import PermissionDenied
from django.db.models import Count, IntegerField, OuterRef, QuerySet, StringAgg, Subquery, TextField, Value
from django.http import HttpRequest
from django.shortcuts import get_object_or_404

from analysis.models import VariantTag
from annotation.models.models import (
    AnnotationVersion,
    GeneAnnotationVersion,
    InvalidAnnotationVersionError,
)
from genes.models import (
    CanonicalTranscript,
    CanonicalTranscriptCollection,
    GeneAnnotationRelease,
    GeneCoverageCanonicalTranscript,
    GeneCoverageCollection,
    GeneList,
    GeneListGeneSymbol,
    GeneSymbol,
    GeneSymbolWiki,
    ReleaseGeneVersion,
    TranscriptVersion,
)
from genes.models_enums import AnnotationConsortium, GeneSymbolAliasSource
from library.utils import pretty_label, update_dict_of_dict_values
from snpdb.grid_columns.custom_columns import get_custom_column_fields_override_and_sample_position
from snpdb.grids import AbstractVariantGrid
from snpdb.models import ImportStatus, Q, Tag, UserSettings, VariantGridColumn
from snpdb.models.models_genome import GenomeBuild
from snpdb.variant_queries import (
    get_variant_queryset_for_gene_symbol,
    variant_qs_filter_has_internal_data,
)
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


class GeneListColumns(DatatableConfig[GeneList]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn('id',
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn('name'),
            RichColumn('user__username', name='user', orderable=True),
            RichColumn('import_status', orderable=True, renderer=self.render_import_status),
            RichColumn('created', client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn('modified', client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
            RichColumn('num_genes'),
        ]

    @staticmethod
    def render_import_status(row: dict[str, Any]):
        return ImportStatus(row["import_status"]).label

    def get_initial_queryset(self) -> QuerySet[GeneList]:
        qs = GeneList.filter_for_user(self.user, success_only=False)
        qs = qs.filter(category__isnull=True)  # only show non-special ones
        qs = qs.select_related("user")  # 'user__username' column - avoid a query per row
        # Count genes via a subquery rather than Count() + GROUP BY: the GROUP BY makes DataTables'
        # unfiltered .count() aggregate every gene list's symbols just to count the lists themselves.
        num_genes_qs = GeneListGeneSymbol.objects.filter(gene_list=OuterRef("pk")) \
            .order_by().values("gene_list").annotate(c=Count("pk")).values("c")
        return qs.annotate(num_genes=Subquery(num_genes_qs, output_field=IntegerField()))

    def filter_queryset(self, qs: QuerySet[GeneList]) -> QuerySet[GeneList]:
        if gene_symbol_id := self.get_query_param('gene_symbol'):  # This is on gene page
            gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol_id)
            qs = GeneList.visible_gene_lists_containing_gene_symbol(qs, gene_symbol)
        return qs


class GeneListGenesColumns(DatatableConfig[GeneListGeneSymbol]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn('original_name', orderable=True),
            RichColumn('gene_symbol', orderable=True, client_renderer="renderGeneSymbol"),
            RichColumn('gene_symbol_alias__source', orderable=True, name='Alias',
                       renderer=self.render_alias_source),
        ]

        self.annotation_releases = {}
        for release in self._get_gene_annotation_releases():
            field_name = f"release_{release.pk}"
            self.annotation_releases[field_name] = release
            self.rich_columns.append(RichColumn(field_name, name=str(release), orderable=True))

    def _get_gene_annotation_releases(self) -> list['GeneAnnotationRelease']:
        return GeneAnnotationRelease.get_for_latest_annotation_versions_for_builds()

    def _get_gene_list(self):
        """ Default is to get list, doing security check against request.user """
        gene_list_id = self.get_query_param("gene_list_id")
        return GeneList.get_for_user(self.user, gene_list_id, success_only=False)

    def get_initial_queryset(self) -> QuerySet[GeneList]:
        gene_list = self._get_gene_list()
        queryset = GeneListGeneSymbol.objects.filter(gene_list=gene_list)

        annotation_kwargs = {}
        for field_name, release in self.annotation_releases.items():
            annotation_kwargs[field_name] = GeneListGeneSymbol.get_joined_genes_qs_annotation_for_release(release)

        return queryset.annotate(**annotation_kwargs)

    @staticmethod
    def render_alias_source(row: dict[str, Any]):
        if alias_source := row["gene_symbol_alias__source"]:
            return GeneSymbolAliasSource(alias_source).label
        return ""


class GeneSymbolVariantsGrid(AbstractVariantGrid):
    """ Uses custom columns subtracting away the gene annotations (as they're displayed above) """
    caption = 'Gene Variants'

    def __init__(self, user, gene_symbol, genome_build_name, **kwargs):
        extra_filters = kwargs.pop("extra_filters", None)
        self.gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol)
        user_settings = UserSettings.get_for_user(user)
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        self.genome_build = genome_build
        self.annotation_version = AnnotationVersion.latest(genome_build)
        fields, override, _ = get_custom_column_fields_override_and_sample_position(user_settings.columns,
                                                                                    self.annotation_version)
        self.fields = self._get_non_gene_fields(fields)
        super().__init__(user)

        af_show_in_percent = settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT
        update_dict_of_dict_values(self._overrides, self._get_standard_overrides(af_show_in_percent))
        update_dict_of_dict_values(self._overrides, override)

        self.extra_filters = extra_filters
        self.extra_config.update({'sortname': "locus__position",
                                  'sortorder': "asc",
                                  'shrinkToFit': False})

    def _get_base_queryset(self) -> QuerySet:
        genes_qs = get_variant_queryset_for_gene_symbol(self.gene_symbol, self.annotation_version)
        return variant_qs_filter_has_internal_data(genes_qs, self.annotation_version, show_clinvar=False)

    def _get_q(self) -> Optional[Q]:
        q_list = []
        if self.extra_filters:
            # Hotspot filters
            if protein_position := self.extra_filters.get("protein_position"):
                transcript_version = TranscriptVersion.objects.get(pk=self.extra_filters["protein_position_transcript_version_id"])
                q_list.append(Q(varianttranscriptannotation__transcript_version=transcript_version,
                                varianttranscriptannotation__protein_position__icontains=protein_position))
            tag_id = self.extra_filters.get("tag")
            if tag_id is not None:  # "" for all tags
                tags_qs = VariantTag.objects.all()
                if tag_id:
                    tag = get_object_or_404(Tag, pk=tag_id)
                    tags_qs = tags_qs.filter(tag=tag)
                alleles_qs = tags_qs.values_list("variant__variantallele__allele")
                q_list.append(Q(variantallele__allele__in=alleles_qs))

        q = None
        if q_list:
            q = reduce(operator.and_, q_list)
        return q

    @staticmethod
    def _get_non_gene_fields(fields):
        """ Remove fields that'll all be the same """
        non_gene_fields = []
        for f in fields:
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


class GenesColumns(DatatableConfig[ReleaseGeneVersion]):
    """ Shows genes / symbols in a gene annotation release """
    grid_name = "Gene Release"
    server_csv_download = True
    search_box_enabled = True
    scroll_x = True
    # The unfiltered count is every gene in the release joined to its annotation, and only feeds the
    # "(filtered from N total)" text
    count_unfiltered = False

    GENE_SYMBOL_FIELD = "gene_version__gene_symbol__symbol"

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        column_labels = _gene_column_labels()
        self.rich_columns = []
        for field in _get_gene_fields():
            is_symbol = field == self.GENE_SYMBOL_FIELD
            self.rich_columns.append(
                RichColumn(key=field, label=column_labels.get(field, pretty_label(field.split("__")[-1])),
                           orderable=True, search=is_symbol,
                           default_sort=SortOrder.ASC if is_symbol else None,
                           client_renderer='renderGeneSymbol' if is_symbol else None))

    def _get_gene_annotation_release_id(self) -> int:
        if release_id := self.get_query_param("gene_annotation_release_id"):
            return int(release_id)
        genome_build = GenomeBuild.get_name_or_alias(self.get_query_param("genome_build_name"))
        av = AnnotationVersion.latest(genome_build)
        return av.gene_annotation_version.gene_annotation_release_id

    def get_initial_queryset(self) -> QuerySet[ReleaseGeneVersion]:
        gene_annotation_release_id = self._get_gene_annotation_release_id()
        gene_annotation_version = GeneAnnotationVersion.objects.filter(
            gene_annotation_release_id=gene_annotation_release_id).order_by("annotation_date").last()
        if gene_annotation_version is None:
            raise InvalidAnnotationVersionError(
                f"No gene annotation version for gene_annotation_release: {gene_annotation_release_id}")

        return ReleaseGeneVersion.objects.filter(
            release_id=gene_annotation_release_id,
            gene_version__gene__geneannotation__version=gene_annotation_version)

    def filter_queryset(self, qs: QuerySet[ReleaseGeneVersion]) -> QuerySet[ReleaseGeneVersion]:
        if contig_id := self.get_query_param("contig_id"):
            gene_versions = TranscriptVersion.objects.filter(contig_id=contig_id).values_list("gene_version_id",
                                                                                              flat=True)
            qs = qs.filter(gene_version__in=gene_versions)

        if column := self.get_query_param("column"):
            if column not in {rc.key for rc in self.enabled_columns}:
                raise PermissionDenied(f"Bad column '{column}'")
            is_null = self.get_query_param("is_null") == "true"
            qs = qs.filter(**{f"{column}__isnull": is_null})
        return qs


def _gene_column_labels() -> dict[str, str]:
    """ Gene column path -> the label the variant grid gives it """
    labels = {}
    for variant_column, label in VariantGridColumn.objects.values_list("variant_column", "label"):
        gene_column = variant_column.replace("variantannotation__", "").replace("transcript_version__", "")
        if gene_column.startswith("gene__"):
            gene_column = "gene_version__" + gene_column
        labels[gene_column] = label
    return labels


class CanonicalTranscriptCollectionColumns(DatatableConfig[CanonicalTranscriptCollection]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn('id',
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn('filename', renderer=self.render_filename),
            RichColumn('description'),
            RichColumn('annotation_consortium', orderable=True, renderer=self.render_annotation_consortium),
            RichColumn('created', client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn('modified', client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
            RichColumn('enrichment_kits'),
            RichColumn('num_transcripts'),
        ]

    @staticmethod
    def render_filename(row: dict[str, Any]):
        filename = row["filename"]
        return os.path.basename(filename)

    @staticmethod
    def render_annotation_consortium(row: dict[str, Any]):
        return AnnotationConsortium(row["annotation_consortium"]).label

    def get_initial_queryset(self) -> QuerySet[GeneList]:
        qs = CanonicalTranscriptCollection.objects.all()
        return qs.annotate(enrichment_kits=StringAgg("enrichmentkit__name", Value(','), output_field=TextField()),
                           num_transcripts=Count("canonicaltranscript"))

    def filter_queryset(self, qs: QuerySet[GeneList]) -> QuerySet[GeneList]:
        if gene_symbol_id := self.get_query_param('gene_symbol'):  # This is on gene page
            gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol_id)
            qs = GeneList.visible_gene_lists_containing_gene_symbol(qs, gene_symbol)
        return qs


class CanonicalTranscriptColumns(DatatableConfig[CanonicalTranscript]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn('gene_symbol__symbol', label="Matched Gene Symbol"),
            RichColumn('transcript__identifier', label="Matched Transcript"),
            RichColumn('original_gene_symbol'),
            RichColumn('original_transcript'),
        ]

    @staticmethod
    def render_annotation_consortium(row: dict[str, Any]):
        return AnnotationConsortium(row["annotation_consortium"]).label

    def get_initial_queryset(self) -> QuerySet[GeneList]:
        return CanonicalTranscript.objects.all()

    def filter_queryset(self, qs: QuerySet[GeneList]) -> QuerySet[GeneList]:
        if collection_id := self.get_query_param('collection_id'):
            qs = qs.filter(collection_id=collection_id)
        return qs


class QCGeneCoverageColumns(DatatableConfig[GeneCoverageCanonicalTranscript]):
    """ Coverage for the genes in a GeneCoverageCollection, optionally restricted to gene lists """
    grid_name = "QC"
    server_csv_download = True

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="gene_symbol__symbol", label="Gene Symbol", orderable=True,
                       client_renderer='renderGeneSymbol'),
            RichColumn(key="transcript__identifier", label="Transcript", orderable=True),
            RichColumn(key="original_gene_symbol", label="Original Symbol", orderable=True),
            RichColumn(key="original_transcript", label="Original Transcript", orderable=True),
            RichColumn(key="min", label="Min", orderable=True, css_class="num"),
            RichColumn(key="mean", label="Mean", orderable=True, css_class="num"),
            RichColumn(key="std_dev", label="Std Dev", orderable=True, css_class="num"),
            RichColumn(key="percent_1x", label="% 1x", orderable=True, css_class="num"),
            RichColumn(key="percent_10x", label="% 10x", orderable=True, css_class="num"),
            RichColumn(key="percent_20x", label="% 20x", orderable=True, css_class="num"),
        ]

    def get_coverage_q(self, gene_coverage_collection, gene_symbols) -> Q:
        filters = [Q(gene_coverage_collection=gene_coverage_collection)]
        if gene_symbols:
            filters.append(Q(gene_symbol__in=gene_symbols))
        return reduce(operator.and_, filters)

    def get_initial_queryset(self) -> QuerySet[GeneCoverageCanonicalTranscript]:
        gene_coverage_collection = get_object_or_404(GeneCoverageCollection,
                                                     pk=self.get_query_param("gene_coverage_collection_id"))
        gene_symbols = set()
        if gene_list_id_list := self.get_query_param("gene_list_id_list"):
            for gene_list_id in gene_list_id_list.split("/"):
                gene_list = GeneList.get_for_user(self.user, gene_list_id, success_only=False)
                gene_symbols.update(gene_list.get_gene_names())

        q = self.get_coverage_q(gene_coverage_collection, gene_symbols)
        return GeneCoverageCanonicalTranscript.objects.filter(q)


class UncoveredGenesColumns(QCGeneCoverageColumns):
    def get_coverage_q(self, gene_coverage_collection, gene_symbols) -> Q:
        q = super().get_coverage_q(gene_coverage_collection, gene_symbols)
        return q & Q(min__lt=int(self.get_query_param("min_depth")))


class GeneSymbolWikiColumns(DatatableConfig[GeneSymbolWiki]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.download_csv_button_enabled = True

        self.rich_columns = [
            RichColumn('gene_symbol', orderable=True, client_renderer="renderGeneSymbol"),
            RichColumn('markdown'),
            RichColumn('last_edited_by__username', name='user', orderable=True),
            RichColumn('created', client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn('modified', client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
        ]

    def get_initial_queryset(self) -> QuerySet[GeneSymbolWiki]:
        return GeneSymbolWiki.objects.all()
