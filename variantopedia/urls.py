from library.django_utils.jqgrid_view import JQGridView
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import perm_path
from variantopedia import views
from variantopedia.grids import AllVariantsGrid, NearbyVariantsGrid, TaggedVariantGrid, \
    VariantTagsGrid, VariantWikiColumns

urlpatterns = [
    perm_path('variants', views.variants, name='variants'),
    perm_path('variants/<genome_build_name>', views.variants, name='genome_build_variants'),

    perm_path('dashboard', views.dashboard, name='dashboard'),
    perm_path('server_status', views.server_status, name='server_status'),
    perm_path('server_status_activity/detail/<int:days_ago>', views.server_status_activity,
              name='server_status_activity_detail'),
    perm_path('server_status_settings/detail', views.server_status_settings, name='server_status_settings_detail'),
    perm_path('database_statistics/detail', views.database_statistics, name='database_statistics_detail'),
    # Tagging
    perm_path('variant_tags/', views.variant_tags, name='variant_tags'),
    perm_path('variant_tags/<genome_build_name>', views.variant_tags, name='genome_build_variant_tags'),

    #
    perm_path('search', views.search, name='search'),
    perm_path('view_variant/<int:variant_id>', views.view_variant, name='view_variant'),
    perm_path('view_variant/<int:variant_id>/<genome_build_name>', views.view_variant,
              name='view_variant_genome_build'),
    perm_path('view_variant_annotation_history/<int:variant_id>', views.view_variant_annotation_history,
              name='view_variant_annotation_history'),
    perm_path('variant_wiki', views.variant_wiki, name='variant_wiki'),
    perm_path('variant_wiki/<genome_build_name>', views.variant_wiki, name='genome_build_variant_wiki'),

    perm_path('view_allele_from_variant/<int:variant_id>', views.view_allele_from_variant,
              name='view_allele_from_variant'),
    perm_path('view_allele/<int:allele_id>', views.view_allele, name='view_allele'),
    perm_path('view_allele/<int:allele_id>/classifications_download', views.export_classifications_allele,
              name='allele_classifications_download'),
    perm_path('allele/<allele_id>/create_variant/<genome_build_name>',
              views.create_variant_for_allele, name='create_variant_for_allele'),
    perm_path('view/<int:variant_id>/<int:annotation_version_id>', views.variant_details_annotation_version,
              name='variant_details_annotation_version'),
    perm_path('nearby_tab/<int:variant_id>/<int:annotation_version_id>', views.nearby_variants_tab,
              name='nearby_variants_tab'),
    perm_path('nearby/<int:variant_id>/<int:annotation_version_id>', views.nearby_variants,
              name='nearby_variants_annotation_version'),
    perm_path('gene_coverage/<slug:gene_symbol_id>', views.gene_coverage, name='gene_coverage'),
    perm_path('variant_sample_information/<int:variant_id>/<genome_build_name>', views.variant_sample_information,
              name='variant_sample_information'),

    # Grids
    perm_path('wiki/datatable/', DatabaseTableView.as_view(column_class=VariantWikiColumns),
              name='variant_wiki_datatable'),
    perm_path('nearby/grid/<variant_id>/<region_type>/<slug:op>/',
              JQGridView.as_view(grid=NearbyVariantsGrid, csv_download=True),
              name='nearby_variants_grid'),
    perm_path('all_variants/grid/<genome_build_name>/<slug:op>/', JQGridView.as_view(grid=AllVariantsGrid, csv_download=True),
              name='all_variants_grid'),
    perm_path('tags/grid/<genome_build_name>/<slug:op>/',
              JQGridView.as_view(grid=VariantTagsGrid, delete_row=True), name='variant_tags_grid'),
    perm_path('tagged_variants/grid/<genome_build_name>/<slug:op>/',
              JQGridView.as_view(grid=TaggedVariantGrid, delete_row=True), name='tagged_variant_grid'),
]
