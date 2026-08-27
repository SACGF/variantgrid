import operator
import re
from functools import cached_property, reduce
from typing import Any, Optional

from django.conf import settings
from django.db import connection
from django.db.models import (
    Case, Count, FilteredRelation, IntegerField, OuterRef, Q, QuerySet, Subquery, TextField, Value, When,
)
from django.db.models.functions import Coalesce, Concat
from django.http import HttpRequest
from django.shortcuts import get_object_or_404
from django.urls import reverse

from analysis.models import Analysis, VariantTag
from annotation.annotation_version_querysets import (
    get_variant_queryset_for_annotation_version,
    get_variant_queryset_for_latest_annotation_version,
)
from annotation.models import AnnotationVersion, VariantAnnotation
from library.utils import JsonDataType, full_class_name, update_dict_of_dict_values
from snpdb.grid_columns.custom_columns import get_custom_column_fields_override_and_sample_position
from snpdb.grids import AbstractVariantGrid, url_if_visible
from snpdb.models import GenomeBuild, Tag, Variant, VariantWiki, VariantZygosityCountCollection
from snpdb.models.models_user_settings import UserGridConfig, UserSettings
from snpdb.utils import get_tag_sort_order_by_tag
from snpdb.variant_filters import get_all_variants_filters, get_variant_filter_q, is_selective
from snpdb.views.datatable_view import CellData, DatatableConfig, RichColumn, SortOrder
from variantopedia.interesting_nearby import get_nearby_qs


def _format_approx_count(n: int) -> str:
    """Format a large approximate count as '~100M', '~1.2B', etc."""
    for threshold, suffix in ((1_000_000_000, 'B'), (1_000_000, 'M'), (1_000, 'K')):
        if n >= threshold:
            rounded = n / threshold
            fmt = f"{rounded:.0f}" if rounded >= 10 else f"{rounded:.1f}"
            return f"~{fmt}{suffix}"
    return f"~{n}"


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
    def render_variant(cell: CellData) -> JsonDataType:
        variant_id = cell["variant"]
        variant = get_object_or_404(Variant, pk=variant_id)
        return {"id": variant_id, "g_hgvs": VariantAnnotation.get_hgvs_g(variant)}

    def render_genome_build(self, _cell: CellData) -> JsonDataType:
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
    # Sorting on a joined or unindexed column full-sorts the whole result set before LIMIT, blowing the
    # statement_timeout (@see issues #1279, #1651). Nothing is user-sortable - every page is served in genomic
    # order (see DEFAULT_ORDER_BY), which is the one ordering a contig-filtered page can stream. @see issue #1663
    SORTABLE_FIELDS: set[str] = set()
    # (contig, position) is the leading edge of the snpdb_locus(contig_id, position, ref_id) unique index, so a
    # contig-filtered page streams straight off it via an incremental sort instead of full-sorting the result set
    # (id-descending walks the whole variant table under a contig filter - measured 100-1000x slower). The pk
    # tiebreaker makes pagination stable. @see issue #1663
    DEFAULT_ORDER_BY = ("locus__contig_id", "locus__position", "pk")

    def __init__(self, user, genome_build_name, **kwargs):
        user_settings = UserSettings.get_for_user(user)
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        self.genome_build = genome_build
        self.annotation_version = AnnotationVersion.latest(genome_build)
        fields, override, _ = get_custom_column_fields_override_and_sample_position(user_settings.columns,
                                                                                    self.annotation_version)
        self.fields = fields
        super().__init__(user)
        af_show_in_percent = settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT
        update_dict_of_dict_values(self._overrides, self._get_standard_overrides(af_show_in_percent))
        update_dict_of_dict_values(self._overrides, override)
        self.vzcc = VariantZygosityCountCollection.get_global_germline_counts()
        self.extra_filters = kwargs.pop("extra_filters", {})
        self.extra_config.update({'sortname': 'locus__position',
                                  'sortorder': "asc",
                                  'shrinkToFit': False})

    def _get_base_queryset(self) -> QuerySet:
        return get_variant_queryset_for_annotation_version(self.annotation_version)

    @cached_property
    def filters(self) -> dict:
        """ The page sends the current selection as extra_filters. A direct grid hit (CSV export, bookmarked
            grid URL) has none, so fall back to what the user last chose on the page """
        if self.extra_filters:
            return self.extra_filters
        return get_all_variants_filters(self.user, self.genome_build)

    def _get_q(self) -> Optional[Q]:
        filters = self.filters
        if not is_selective(filters):
            # Match nothing rather than scanning the whole variant table - the page explains why
            return Q(pk__isnull=True)

        filter_list = [
            get_variant_filter_q(self.genome_build, self.annotation_version,
                                 contig_ids=filters.get("contig_ids"),
                                 non_standard_contigs=filters.get("non_standard_contigs", False),
                                 gene_symbols=filters.get("gene_symbols"),
                                 variant_types=filters.get("variant_types")),
            Variant.get_no_reference_q(),
        ]
        if min_count := int(filters.get("min_count") or 0):
            # benchmarking - I found it much faster to do both of these queries (seems redundant)
            hom_nonzero = Q(**{f"{self.vzcc.hom_alias}__gt": 0})
            het_nonzero = Q(**{f"{self.vzcc.het_alias}__gt": 0})
            filter_list.append(hom_nonzero | het_nonzero)
            filter_list.append(Q(**{f"{self.vzcc.non_ref_call_alias}__gte": min_count}))

        return reduce(operator.and_, filter_list)

    def get_colmodels(self, remove_server_side_only=False):
        """ Only the allowlisted columns keep their sort arrows """
        colmodels = super().get_colmodels(remove_server_side_only=remove_server_side_only)
        for cm in colmodels:
            if cm.get("name") not in self.SORTABLE_FIELDS:
                cm["sortable"] = False
        return colmodels

    def _sort_items(self, items, sidx, sord):
        """ Serve every page in genomic order regardless of any sidx a hand-crafted grid URL supplies. Emitted as
            a plain order_by so it matches the snpdb_locus(contig_id, position, ref_id) btree exactly - the base
            class's F(sidx).asc(nulls_first=...) path defeats that index. @see issue #1663 """
        return items.order_by(*self.DEFAULT_ORDER_BY)

    def _get_approx_count(self, qs) -> int:
        sql, params = qs.query.sql_with_params()
        with connection.cursor() as cursor:
            cursor.execute(f"EXPLAIN {sql}", params)
            first_line = cursor.fetchone()[0]
        match = re.search(r'rows=(\d+)', first_line)
        if not match:
            raise ValueError(f"Could not parse row estimate from EXPLAIN output: {first_line!r}")
        return int(match.group(1))

    def get_known_count(self, request, items) -> Optional[int]:
        """ A COUNT(*) over a huge table costs more than the page itself - hand the paginator the
            planner's estimate instead, and tell the user it's approximate """
        if self.get_filters(request):
            return None  # jqGrid column filters narrow the rows the estimate was taken over

        try:
            estimate = self._get_approx_count(items)
        except Exception:
            return None

        if estimate >= 1_000_000:
            self._used_approx_count = True
            return estimate
        return None

    def get_data(self, request) -> dict:
        self._used_approx_count = False
        data = super().get_data(request)
        if self._used_approx_count:
            data['approximate_records'] = _format_approx_count(data['records'])
        return data


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


class VariantTagsColumns(DatatableConfig[VariantTag]):
    """ List VariantTags (Tag-centric) - @see variant_tags.html and base_related_analyses.html """
    GRID_NAME = 'Variant Tags'

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.genome_build = GenomeBuild.get_name_or_alias(self.get_query_param("genome_build_name"))

        self.rich_columns = [
            RichColumn("id", visible=False),
            RichColumn("variant_string", label="Variant", orderable=True, default_sort=SortOrder.ASC,
                       extra_columns=["id", "variant__id", "tag__id", "analysis__id"],
                       renderer=self.render_variant, client_renderer="renderVariantTagVariant"),
            RichColumn(name="genome_build", label="Genome Build", renderer=self.render_genome_build),
            RichColumn("variant__variantannotation__transcript_version__gene_version__gene_symbol__symbol",
                       name="gene_symbol", label="Gene", orderable=True,
                       client_renderer="renderGeneSymbolNewWindow"),
            RichColumn("tag__id", name="tag", label="Tag", orderable=True, extra_columns=["variant__id"],
                       renderer=self.render_tag, client_renderer="renderVariantTagPill"),
            RichColumn("analysis__name", name="analysis", label="Analysis", orderable=True,
                       extra_columns=["analysis__id"],
                       renderer=self.render_analysis, client_renderer="renderVariantTagAnalysis"),
            RichColumn("user__username", name="user", label="Username", orderable=True),
            RichColumn("created", label="Created", orderable=True, client_renderer="TableFormat.timestamp"),
            RichColumn("id", name="delete", label="", extra_columns=["analysis__id"],
                       renderer=self.render_delete, client_renderer="TableFormat.deleteRow"),
        ]

    def render_genome_build(self, _cell: CellData) -> JsonDataType:
        return self.genome_build.name

    @staticmethod
    def render_variant(cell: CellData) -> JsonDataType:
        data = {
            "variant_string": cell.value,
            "url": url_if_visible("view_variant", variant_id=cell["variant__id"]),
        }
        # The tag is a to-do item - offer to complete it @see analysis.variant_tag_operations
        if cell["tag__id"] == settings.TAG_REQUIRES_CLASSIFICATION and (analysis_id := cell["analysis__id"]):
            data["classify_url"] = url_if_visible("create_classification_for_variant_tag",
                                                  analysis_id=analysis_id, variant_tag_id=cell["id"])
        return data

    @staticmethod
    def render_tag(cell: CellData) -> JsonDataType:
        return {"tag": cell.value, "variant_id": cell["variant__id"]}

    @staticmethod
    def render_analysis(cell: CellData) -> JsonDataType:
        analysis_id = cell["analysis__id"]
        if analysis_id is None:
            return None
        return {
            "text": f"{analysis_id} - {cell.value}",
            "url": url_if_visible("analysis", analysis_id=analysis_id),
        }

    def render_delete(self, cell: CellData) -> Optional[str]:
        """ can_write delegates to the analysis, so build the stub from the row rather than re-loading the tag """
        variant_tag = VariantTag(pk=cell.value, analysis_id=cell["analysis__id"])
        if not variant_tag.can_write(self.user):
            return None
        return reverse('group_permissions_object_delete',
                       kwargs={'class_name': full_class_name(VariantTag), 'primary_key': cell.value})

    def get_initial_queryset(self) -> QuerySet[VariantTag]:
        # get_for_build has already restricted this to tags visible in the build, either via their
        # allele or - for a tag made in this build - via the tag's own variant
        qs = VariantTag.get_for_build(self.genome_build)
        # Pick the variant for *this* build out of the allele, rather than whichever one a plain join
        # would land on
        qs = qs.annotate(build_variant_allele=FilteredRelation(
            "allele__variantallele",
            condition=Q(allele__variantallele__genome_build=self.genome_build)))
        return self._annotate_variant_string(qs)

    @staticmethod
    def _annotate_variant_string(qs: QuerySet[VariantTag]) -> QuerySet[VariantTag]:
        """ A "1:123321 G>C" string built from the build's variant, falling back per field to the one
            the tag was made on - a tag keeps its own variant until the liftover task assigns it an
            allele (@see analysis.tasks.variant_tag_tasks._liftover_variant_tag), and that variant is
            in this build by definition. Dropping those rows made the grid disagree with the tag counts.

            The fallback is per field because Concat renders a NULL argument as empty rather than
            returning NULL, so coalescing the finished strings would never reach the second one. """
        def field(name: str):
            return Coalesce(f"build_variant_allele__variant__{name}", f"variant__{name}")

        return qs.annotate(variant_string=Concat(
            field("locus__contig__name"), Value(":"), field("locus__position"), Value(" "),
            field("locus__ref__seq"), Value(">"), field("alt__seq"), output_field=TextField()))

    def filter_queryset(self, qs: QuerySet[VariantTag]) -> QuerySet[VariantTag]:
        analysis_ids = self.get_query_json("analysis_ids")
        if analysis_ids is not None:
            analyses_queryset = Analysis.filter_for_user(self.user).filter(pk__in=analysis_ids)
            qs = qs.filter(analysis__in=analyses_queryset)

        if gene_id := self.get_query_param("gene"):
            qs = qs.filter(variant__variantannotation__transcript_version__gene_version__gene_id=gene_id)

        if tag_id := self.get_query_param("tag"):
            qs = qs.filter(tag_id=tag_id)

        if tag_ids := self.get_query_json("tags"):
            qs = qs.filter(tag__in=tag_ids)

        filter_user_id = self.get_query_param("user")
        if filter_user_id:
            qs = qs.filter(user_id=filter_user_id)

        user_grid_config = UserGridConfig.get(self.user, self.GRID_NAME)
        if user_grid_config.show_group_data or filter_user_id:
            # An explicit user filter overrides show_group_data - still permission checked
            qs = VariantTag.filter_for_user(self.user, queryset=qs)
        else:
            qs = qs.filter(user=self.user)
        return qs


class TaggedVariantGrid(AbstractVariantGrid):
    """ Shows Variants that have been tagged (Variant-centric) """
    caption = 'Variant with tags'

    TAG_COUNT_OVERRIDE = {
        'model_field': False, 'queryset_field': False,
        'name': 'tag_count', 'index': 'tag_count', 'label': 'Tag Events',
        'width': 60, 'sorttype': 'int',
    }

    def __init__(self, user, genome_build_name, extra_filters=None):
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        self.genome_build = genome_build
        tag_ids = []
        require_all_tags = False
        filter_user_id = None
        if extra_filters:
            if tag_id := extra_filters.get("tag"):
                tag_ids.append(tag_id)
            # Variants carrying ALL of these tags - @see tag stats co-occurrence card
            if all_tag_ids := extra_filters.get("tags"):
                tag_ids.extend(all_tag_ids)
                require_all_tags = True
            filter_user_id = extra_filters.get("user")
        self.tag_ids = tag_ids
        self.require_all_tags = require_all_tags
        self.filter_user_id = filter_user_id

        user_settings = UserSettings.get_for_user(user)
        self.annotation_version = AnnotationVersion.latest(genome_build)
        fields, override, _ = get_custom_column_fields_override_and_sample_position(user_settings.columns,
                                                                                    self.annotation_version)
        self.fields = fields + ["tag_count"]
        super().__init__(user)

        af_show_in_percent = settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT
        update_dict_of_dict_values(self._overrides, self._get_standard_overrides(af_show_in_percent))
        update_dict_of_dict_values(self._overrides, override)
        update_dict_of_dict_values(self._overrides, {"tag_count": self.TAG_COUNT_OVERRIDE})
        self.extra_config.update({'sortname': "locus__position",
                                  'sortorder': "asc",
                                  'shrinkToFit': False})

    def _get_grid_only_annotation_kwargs(self):
        """ How many times this variant has been tagged - sort on it to find the most re-tagged variants """
        a_kwargs = super()._get_grid_only_annotation_kwargs()
        tag_count_qs = VariantTag.filter_for_user(self.user).filter(
            allele__variantallele__variant_id=OuterRef("id")).values("allele").annotate(
            tag_count=Count("pk")).values_list("tag_count")
        a_kwargs["tag_count"] = Subquery(tag_count_qs[:1])
        return a_kwargs

    def _get_base_queryset(self) -> QuerySet:
        genome_build = self.annotation_version.genome_build
        qs = get_variant_queryset_for_latest_annotation_version(genome_build)
        qs = qs.filter(Variant.get_contigs_q(genome_build))
        return qs

    def _get_q(self) -> Optional[Q]:
        genome_build = self.annotation_version.genome_build
        user_grid_config = UserGridConfig.get(self.user, self.caption)
        tags_qs = VariantTag.filter_for_user(self.user)
        if self.filter_user_id:
            # An explicit user filter overrides show_group_data - still permission checked
            tags_qs = tags_qs.filter(user_id=self.filter_user_id)
        elif not user_grid_config.show_group_data:
            tags_qs = tags_qs.filter(user=self.user)
        if self.require_all_tags:
            q_list = [VariantTag.variants_for_build_q(genome_build, tags_qs, [tag_id]) for tag_id in self.tag_ids]
            return reduce(operator.and_, q_list)
        return VariantTag.variants_for_build_q(genome_build, tags_qs, self.tag_ids)


class VariantTagCountsColumns(DatatableConfig[VariantTag]):
    """ This is for showing on the variant page """
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.expand_client_renderer = DatatableConfig._row_expand_ajax('viewTagDetail', id_field="tag")
        self.tag_stale_date = UserSettings.get_for_user(self.user).variant_tag_stale_date
        # Custom tag ordering from the tag colors page - tags without an entry sort as 0
        self.sort_order_by_tag = get_tag_sort_order_by_tag(self.user)
        tag_sort_keys = ['tag_sort_order', 'tag'] if self.sort_order_by_tag else ['tag']

        self.rich_columns = [
            RichColumn('tag', client_renderer='tagRenderer', sort_keys=tag_sort_keys,
                       default_sort=SortOrder.ASC),
            RichColumn('count', orderable=True),
        ]
        if self.tag_stale_date:
            self.rich_columns.append(RichColumn('fresh_count', label='Fresh', orderable=True))
        self.rich_columns += [
            RichColumn('last_created', client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn('last_created', name='time_ago', client_renderer='TableFormat.timeAgo'),
        ]

    def _get_sort_tiebreaker(self) -> str:
        return "tag"

    def get_initial_queryset(self) -> QuerySet[VariantTag]:
        variant_id = self.get_query_param('variant_id')
        variant = Variant.objects.get(pk=variant_id)
        qs = VariantTag.get_variant_tag_counts_qs(variant)
        if self.tag_stale_date:
            qs = qs.annotate(fresh_count=Count("id", filter=Q(created__gte=self.tag_stale_date)))
        if self.sort_order_by_tag:
            whens = [When(tag=tag_id, then=Value(sort_order))
                     for tag_id, sort_order in self.sort_order_by_tag.items()]
            qs = qs.annotate(tag_sort_order=Case(*whens, default=Value(0), output_field=IntegerField()))
        return qs


class VariantTagDetailColumns(DatatableConfig[VariantTag]):
    """ This is the detail expanded on variant tags page """
    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.rich_columns = [
            RichColumn('id', client_renderer='tagDetailRenderer'),
            RichColumn('id', name='can_write', visible=False, renderer=self.can_write),
            RichColumn('analysis', client_renderer='analysisLinkRenderer'),
            RichColumn('user__username', name='user', orderable=True),
            RichColumn('created', client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
            RichColumn('created', name='time_ago', client_renderer='TableFormat.timeAgo'),
        ]

    def can_write(self, row: dict[str, Any]) -> bool:
        """ This is really inefficient as it instantiates an object per row that has already had values() called on
            it. Perhaps it would be more efficient to be able to swap in a serializer to produce each row
            however that takes a lot of messing around with DatabaseTableView/DatatableConfig """
        variant_tag = VariantTag.objects.get(pk=row["id"])
        return variant_tag.can_write(self.user)

    def get_initial_queryset(self) -> QuerySet[VariantTag]:
        variant_id = self.get_query_param('variant_id')
        tag_name = self.get_query_param('tag')

        variant = Variant.objects.get(pk=variant_id)
        tag = Tag.objects.get(pk=tag_name)
        # Not going to use anything build specific so don't care about build
        genome_build = variant.any_genome_build
        qs = VariantTag.get_for_build(genome_build, variant_qs=variant.equivalent_variants)
        qs = VariantTag.filter_for_user(self.user, queryset=qs)
        qs = qs.filter(tag=tag)
        return qs
