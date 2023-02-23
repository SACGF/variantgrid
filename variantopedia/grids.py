import operator
from functools import reduce
from typing import Dict, Any, Optional

from django.conf import settings
from django.db.models import TextField, Value, QuerySet, Q
from django.http import HttpRequest
from django.shortcuts import get_object_or_404

from analysis.models import VariantTag, Analysis
from annotation.annotation_version_querysets import get_variant_queryset_for_latest_annotation_version, \
    get_variant_queryset_for_annotation_version
from annotation.models import AnnotationVersion
from genes.hgvs import HGVSMatcher
from library.jqgrid.jqgrid_user_row_config import JqGridUserRowConfig
from library.utils import update_dict_of_dict_values
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

    def render_genome_build(self, _row: Dict[str, Any]) -> JsonDataType:
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

    def __init__(self, user, genome_build_name, **kwargs):
        user_settings = UserSettings.get_for_user(user)
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        self.annotation_version = AnnotationVersion.latest(genome_build)
        fields, override, _ = get_custom_column_fields_override_and_sample_position(user_settings.columns,
                                                                                    self.annotation_version)
        fields.remove("tags")
        self.fields = fields
        super().__init__(user)
        af_show_in_percent = settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT
        update_dict_of_dict_values(self._overrides, self._get_standard_overrides(af_show_in_percent))
        update_dict_of_dict_values(self._overrides, override)
        self.vzcc = VariantZygosityCountCollection.get_global_germline_counts()
        self.extra_filters = kwargs.pop("extra_filters", {})
        self.extra_config.update({'sortname': self.vzcc.germline_counts_alias,
                                  'sortorder': "desc",
                                  'shrinkToFit': False})

    def _get_base_queryset(self) -> QuerySet:
        return get_variant_queryset_for_annotation_version(self.annotation_version)

    def _get_q(self) -> Optional[Q]:
        filter_list = [
            Variant.get_contigs_q(self.annotation_version.genome_build),
        ]
        if self.extra_filters:
            if min_count := int(self.extra_filters.get("min_count", 0)):
                filter_list.append(Q(**{f"{self.vzcc.germline_counts_alias}__gte": min_count}))
        else:
            # benchmarking - I found it much faster to do both of these queries (seems redundant)
            hom_nonzero = Q(**{f"{self.vzcc.hom_alias}__gt": 0})
            het_nonzero = Q(**{f"{self.vzcc.het_alias}__gt": 0})
            filter_list.append(hom_nonzero | het_nonzero)
            filter_list.append(Q(**{f"{self.vzcc.germline_counts_alias}__gt": 0}))

        return reduce(operator.and_, filter_list)


class NearbyVariantsGrid(AbstractVariantGrid):
    caption = 'Nearby Variants'

    def __init__(self, user, variant_id, genome_build_name, region_type, gene_symbol=None, **kwargs):
        self.variant = get_object_or_404(Variant, pk=variant_id)
        self.genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        self.region_type = region_type
        self.gene_symbol = gene_symbol

        user_settings = UserSettings.get_for_user(user)
        self.annotation_version = AnnotationVersion.latest(self.genome_build)
        fields, override, _ = get_custom_column_fields_override_and_sample_position(user_settings.columns,
                                                                                    self.annotation_version)
        fields.remove("tags")
        self.fields = fields
        super().__init__(user)
        af_show_in_percent = settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT
        update_dict_of_dict_values(self._overrides, self._get_standard_overrides(af_show_in_percent))
        update_dict_of_dict_values(self._overrides, override)
        self.extra_config.update({'sortname': "locus__position",
                                  'sortorder': "desc",
                                  'shrinkToFit': False})

    def _get_base_queryset(self) -> QuerySet:
        region_filters = get_nearby_qs(self.variant, self.annotation_version)
        rf_data = region_filters[self.region_type]
        if self.gene_symbol:
            qs = rf_data[self.gene_symbol]
        else:
            qs = rf_data
        return qs


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
            {'index': 'view_genome_build', 'name': 'view_genome_build', 'label': 'Genome Build', 'sortable': False},
        ]
        colmodels = super().get_colmodels(remove_server_side_only=remove_server_side_only)
        return before_colmodels + colmodels


class TaggedVariantGrid(AbstractVariantGrid):
    """ Shows Variants that have been tagged (Variant-centric) """
    caption = 'Variant with tags'

    def __init__(self, user, genome_build_name, extra_filters=None):
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        tag_ids = []
        if extra_filters:
            if tag_id := extra_filters.get("tag"):
                tag_ids.append(tag_id)
        self.tag_ids = tag_ids

        user_settings = UserSettings.get_for_user(user)
        self.annotation_version = AnnotationVersion.latest(genome_build)
        fields, override, _ = get_custom_column_fields_override_and_sample_position(user_settings.columns,
                                                                                    self.annotation_version)
        fields.remove("tags")
        self.fields = fields
        super().__init__(user)

        af_show_in_percent = settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT
        update_dict_of_dict_values(self._overrides, self._get_standard_overrides(af_show_in_percent))
        update_dict_of_dict_values(self._overrides, override)
        self.extra_config.update({'sortname': "locus__position",
                                  'sortorder': "asc",
                                  'shrinkToFit': False})

    def _get_base_queryset(self) -> QuerySet:
        genome_build = self.annotation_version.genome_build
        qs = get_variant_queryset_for_latest_annotation_version(genome_build)
        qs = qs.filter(Variant.get_contigs_q(genome_build))
        return qs

    def _get_q(self) -> Optional[Q]:
        genome_build = self.annotation_version.genome_build
        user_grid_config = UserGridConfig.get(self.user, self.caption)
        tags_qs = VariantTag.filter_for_user(self.user)
        if not user_grid_config.show_group_data:
            tags_qs = tags_qs.filter(user=self.user)
        return VariantTag.variants_for_build_q(genome_build, tags_qs, self.tag_ids)
