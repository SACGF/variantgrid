from typing import Dict, Any

from django.db.models import TextField, Value, QuerySet
from django.http import HttpRequest
from django.shortcuts import get_object_or_404

from analysis.models import VariantTag, Analysis
from annotation.annotation_version_querysets import get_variant_queryset_for_latest_annotation_version
from annotation.models import AnnotationVersion
from genes.hgvs import HGVSMatcher
from library.jqgrid_user_row_config import JqGridUserRowConfig
from snpdb.grid_columns.custom_columns import get_custom_column_fields_override_and_sample_position
from snpdb.grids import AbstractVariantGrid
from snpdb.models import Variant, VariantZygosityCountCollection, GenomeBuild, Tag, VariantWiki
from snpdb.models.models_user_settings import UserSettings, UserGridConfig
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder
from uicore.json.json_types import JsonDataType
from variantopedia.interesting_nearby import get_nearby_qs


class VariantWikiColumns(DatatableConfig[VariantWiki]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.download_csv_button_enabled = True

        # self.expand_client_renderer = DatatableConfig._row_expand_ajax('eventlog_detail', expected_height=120)
        self.rich_columns = [
            RichColumn('variant', renderer=self.render_variant, client_renderer="renderVariantId"),
            RichColumn(name="genome_build", renderer=self.render_genome_build, visible=False),
            RichColumn('markdown'),
            RichColumn('last_edited_by__username', name='user', orderable=True),
            RichColumn('created', client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn('modified', client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
        ]

    @staticmethod
    def render_variant(row: Dict[str, Any]) -> JsonDataType:
        variant_id = row["variant"]
        variant = get_object_or_404(Variant, pk=variant_id)
        genome_build = next(iter(variant.genome_builds))
        g_hgvs = HGVSMatcher(genome_build).variant_to_g_hgvs(variant)
        return {"id": variant_id, "g_hgvs": g_hgvs}

    def render_genome_build(self, row: Dict[str, Any]) -> JsonDataType:
        return self.get_query_param('genome_build')

    def get_initial_queryset(self) -> QuerySet[VariantWiki]:
        return VariantWiki.objects.all()

    def filter_queryset(self, qs: QuerySet[VariantWiki]) -> QuerySet[VariantWiki]:
        if genome_build_name := self.get_query_param('genome_build'):
            genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
            qs = qs.filter(variant__locus__contig__genomebuildcontig__genome_build=genome_build)
        return qs


class AllVariantsGrid(AbstractVariantGrid):
    caption = 'All Variants'
    fields = ["id", "locus__contig__name", 'locus__position', 'locus__ref', 'alt']
    colmodel_overrides = {
        'id': {'editable': False, 'width': 90, 'fixed': True, 'formatter': 'detailsLink'},
        'tags_global': {'classes': 'no-word-wrap', 'formatter': 'tagsGlobalFormatter', 'sortable': False},
    }

    def __init__(self, user, **kwargs):
        user_settings = UserSettings.get_for_user(user)
        fields, override, _ = get_custom_column_fields_override_and_sample_position(user_settings.columns)
        self.fields = fields
        super().__init__(user)
        self.update_overrides(override)

        extra_filters = kwargs.pop("extra_filters", {})
        queryset = get_variant_queryset_for_latest_annotation_version(user_settings.default_genome_build)
        queryset, count_column = VariantZygosityCountCollection.annotate_global_germline_counts(queryset)

        filter_kwargs = {}
        if extra_filters:
            try:
                min_count = int(extra_filters.get("min_count"))
            except TypeError:
                min_count = 0
            filter_kwargs[count_column + "__gte"] = min_count
        else:
            filter_kwargs[count_column + "__gt"] = 0  # By default only show those that have samples

        queryset = queryset.filter(**filter_kwargs)
        self.queryset = queryset.values(*self.get_queryset_field_names())
        self.extra_config.update({'sortname': count_column,
                                  'sortorder': "desc",
                                  'shrinkToFit': False})


class NearbyVariantsGrid(AbstractVariantGrid):
    caption = 'Nearby Variants'
    fields = ["id", "locus__contig__name", 'locus__position', 'locus__ref', 'alt']
    colmodel_overrides = {'id': {'editable': False, 'width': 90, 'fixed': True, 'formatter': 'detailsLink'}}

    def __init__(self, user, variant_id, region_type, **kwargs):
        super().__init__(user)

        variant = get_object_or_404(Variant, pk=variant_id)

        user_settings = UserSettings.get_for_user(user)
        fields, override, _ = get_custom_column_fields_override_and_sample_position(user_settings.columns)
        self.fields = fields
        self.update_overrides(override)

        annotation_version = AnnotationVersion.latest(user_settings.default_genome_build)
        region_filters = get_nearby_qs(variant, annotation_version)
        queryset = region_filters[region_type]
        self.queryset = queryset.values(*self.get_queryset_field_names())
        self.extra_config.update({'sortname': "locus__position",
                                  'sortorder': "desc",
                                  'shrinkToFit': False})


class VariantTagsGrid(JqGridUserRowConfig):
    """ List VariantTags (Tag-centric) """
    model = VariantTag
    caption = 'Variant Tags'
    fields = ["id", "variant__variantannotation__transcript_version__gene_version__gene_symbol__symbol",
              "variant__id", "node__id", "tag__id", "analysis__name", "analysis__id", "user__username", "created"]

    colmodel_overrides = {
        'id': {'hidden': True},
        "variant__id": {"hidden": True},
        "node__id": {"hidden": True},
        "variant__variantannotation__transcript_version__gene_version__gene_symbol__symbol": {'label': 'Gene', 'formatter': 'geneSymbolNewWindowLink'},
        "tag__id": {'label': "Tag", "formatter": "formatVariantTag"},
        "analysis__name": {'label': 'Analysis', "formatter": "formatAnalysis"},
        "analysis__id": {'hidden': True},
        "user__username": {'label': "Username"},
        "created": {'label': "Created"},
    }

    def __init__(self, user, genome_build_name, extra_filters=None, **kwargs):
        super().__init__(user)

        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        queryset = VariantTag.get_for_build(genome_build)

        if extra_filters:
            analysis_ids = extra_filters.get("analysis_ids")
            if analysis_ids is not None:
                analyses_queryset = Analysis.filter_for_user(user).filter(pk__in=analysis_ids)
                queryset = queryset.filter(analysis__in=analyses_queryset)

            gene_id = extra_filters.get("gene")
            if gene_id:
                queryset = queryset.filter(variant__variantannotation__transcript_version__gene_version__gene_id=gene_id)

            tag_id = extra_filters.get("tag")
            if tag_id is not None:
                tag = Tag.objects.get(pk=tag_id)
                queryset = queryset.filter(tag=tag)

        user_grid_config = UserGridConfig.get(user, self.caption)
        if user_grid_config.show_group_data:
            queryset = VariantTag.filter_for_user(user, queryset=queryset)
        else:
            queryset = queryset.filter(user=user)

        # Need to go through Allele to get variant in this build
        queryset = queryset.filter(allele__variantallele__genome_build=genome_build)
        queryset = Variant.annotate_variant_string(queryset,
                                                   path_to_variant="allele__variantallele__variant__")
        queryset = queryset.annotate(view_genome_build=Value(genome_build_name, output_field=TextField()))
        field_names = self.get_field_names() + ["variant_string", "view_genome_build"]
        self.queryset = queryset.values(*field_names)
        self.extra_config.update({'sortname': 'variant_string',
                                  'sortorder': 'asc'})

    def get_colmodels(self, remove_server_side_only=False):
        before_colmodels = [
            {'index': 'variant_string', 'name': 'variant_string',
             'label': 'Variant', 'formatter': 'formatVariantTagFirstColumn'},
            {'index': 'view_genome_build', 'name': 'view_genome_build', 'label': 'Genome Build'},
        ]
        colmodels = super().get_colmodels(remove_server_side_only=remove_server_side_only)
        return before_colmodels + colmodels


class TaggedVariantGrid(AbstractVariantGrid):
    """ Shows Variants that have been tagged (Variant-centric) """
    model = Variant
    caption = 'Variant with tags'
    fields = ["id", "locus__contig__name", 'locus__position', 'locus__ref', 'alt']
    colmodel_overrides = {
        'id': {'editable': False, 'width': 90, 'fixed': True, 'formatter': 'detailsLink'},
        'tags_global': {'classes': 'no-word-wrap', 'formatter': 'tagsGlobalFormatter', 'sortable': False},
        "variantannotation__transcript_version__gene_version__gene_symbol__symbol": {'label': 'Gene',
                                                                                     'formatter': 'geneSymbolNewWindowLink'},
    }

    def __init__(self, user, genome_build_name, extra_filters=None):
        super().__init__(user)

        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        user_grid_config = UserGridConfig.get(user, self.caption)
        tags_qs = VariantTag.filter_for_user(user)
        if not user_grid_config.show_group_data:
            tags_qs = tags_qs.filter(user=user)

        tag_ids = []
        if extra_filters:
            if tag_id := extra_filters.get("tag"):
                tag_ids.append(tag_id)

        qs = get_variant_queryset_for_latest_annotation_version(genome_build)
        qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
        qs = qs.filter(VariantTag.variants_for_build_q(genome_build, tags_qs, tag_ids))
        user_settings = UserSettings.get_for_user(user)
        fields, override, _ = get_custom_column_fields_override_and_sample_position(user_settings.columns)
        fields.remove("tags")
        self.fields = fields
        self.update_overrides(override)

        self.queryset = qs.values(*self.get_queryset_field_names())
        self.extra_config.update({'sortname': "locus__position",
                                  'sortorder': "asc",
                                  'shrinkToFit': False})
