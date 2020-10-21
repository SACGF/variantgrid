from django.shortcuts import get_object_or_404
from django.urls.base import reverse

from annotation.annotation_version_querysets import get_variant_queryset_for_latest_annotation_version
from library.django_utils.jqgrid_view import JQGridViewOp
from library.jqgrid_user_row_config import JqGridUserRowConfig
from snpdb.grid_columns.custom_columns import get_custom_column_fields_override_and_sample_position
from snpdb.grids import AbstractVariantGrid
from snpdb.models import Variant, VariantZygosityCountCollection
from snpdb.models.models_user_settings import UserSettings
from variantopedia.interesting_nearby import filter_variant_exon, filter_variant_codon, filter_variant_domain, \
    filter_variant_range


class VariantWikiGrid(JqGridUserRowConfig):
    model = Variant
    caption = 'Variant Wiki'
    fields = ["id", 'variantwiki__markdown', 'variantwiki__last_edited_by__username']
    colmodel_overrides = {'id': {'width': 20, 'formatter': 'viewVariantDetails'},
                          'variantwiki__markdown': {'label': 'Markdown', 'formatter': 'viewContent'},
                          'variantwiki__last_edited_by__username': {'label': 'Last edited by'}}

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.filter(variantwiki__isnull=False)
        self.queryset = queryset.values(*self.get_field_names())
        grid_export_url = reverse("variantopedia_wiki_grid", kwargs={"op": JQGridViewOp.DOWNLOAD})
        self.extra_config['grid_export_url'] = grid_export_url


class AllVariantsGrid(AbstractVariantGrid):
    caption = 'All Variants'
    fields = ["id", "locus__contig__name", 'locus__position', 'locus__ref', 'alt']
    colmodel_overrides = {'id': {'editable': False, 'width': 90, 'fixed': True, 'formatter': 'detailsLink'}}

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
        variant = get_object_or_404(Variant, pk=variant_id)

        user_settings = UserSettings.get_for_user(user)
        fields, override, _ = get_custom_column_fields_override_and_sample_position(user_settings.columns)
        self.fields = fields
        super().__init__(user)
        self.update_overrides(override)

        queryset = get_variant_queryset_for_latest_annotation_version(user_settings.default_genome_build)
        queryset, count_column = VariantZygosityCountCollection.annotate_global_germline_counts(queryset)

        REGION_FILTERS = {
            "codon": filter_variant_codon,
            "exon": filter_variant_exon,
            "domain": filter_variant_domain,
            "range": filter_variant_range,
        }
        filter_func = REGION_FILTERS[region_type]
        queryset = filter_func(queryset, variant)
        self.queryset = queryset.values(*self.get_queryset_field_names())
        self.extra_config.update({'sortname': count_column,
                                  'sortorder': "desc",
                                  'shrinkToFit': False})
