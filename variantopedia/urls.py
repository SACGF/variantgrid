from library.django_utils.jqgrid_view import JQGridView
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path
from variantopedia import views
from variantopedia.grids import AllVariantsGrid, NearbyVariantsGrid, TaggedVariantGrid, \
    VariantTagsGrid, VariantWikiColumns, VariantTagCountsColumns, VariantTagDetailColumns

urlpatterns = [
    path('variants', views.variants, name='variants'),
    path('variants/<genome_build_name>', views.variants, name='genome_build_variants'),

    path('dashboard', views.dashboard, name='dashboard'),
    path('server_status', views.server_status, name='server_status'),
    path('server_status_activity/detail/<int:days_ago>', views.server_status_activity,
         name='server_status_activity_detail'),
    path('server_status_settings/detail', views.server_status_settings, name='server_status_settings_detail'),
    path('health_check_details', views.health_check_details, name='health_check_details'),
    path('database_statistics/detail', views.database_statistics, name='database_statistics_detail'),
    # Tagging
    path('variant_tags/', views.variant_tags, name='variant_tags'),
    path('variant_tags/<genome_build_name>', views.variant_tags, name='genome_build_variant_tags'),

    #
    path('search', views.search, name='search'),
    path('view_variant/<int:variant_id>', views.view_variant, name='view_variant'),
    path('view_variant/<int:variant_id>/<genome_build_name>', views.view_variant,
         name='view_variant_genome_build'),
    path('view_variant_annotation_history/<int:variant_id>', views.view_variant_annotation_history,
         name='view_variant_annotation_history'),
    path('variant_wiki', views.variant_wiki, name='variant_wiki'),
    path('variant_wiki/<genome_build_name>', views.variant_wiki, name='genome_build_variant_wiki'),

    path('view_allele_from_variant/<int:variant_id>', views.view_allele_from_variant,
         name='view_allele_from_variant'),
    path('view_allele/<int:allele_id>', views.view_allele, name='view_allele'),
    path('a<int:allele_id>', views.view_allele, name='view_allele_compact'),
    path('view_allele/<int:allele_id>/classifications_download', views.export_classifications_allele,
         name='allele_classifications_download'),
    path('allele/<allele_id>/create_variant/<genome_build_name>',
         views.create_variant_for_allele, name='create_variant_for_allele'),
    path('view/<int:variant_id>/<int:annotation_version_id>', views.variant_details_annotation_version,
         name='variant_details_annotation_version'),
    path('nearby_tab/<int:variant_id>/<int:annotation_version_id>', views.nearby_variants_tab,
         name='nearby_variants_tab'),
    path('nearby/<int:variant_id>/<int:annotation_version_id>', views.nearby_variants,
         name='nearby_variants_annotation_version'),
    path('gene_coverage/<slug:gene_symbol_id>', views.gene_coverage, name='gene_coverage'),
    path('variant_sample_information/<int:variant_id>/<genome_build_name>', views.variant_sample_information,
         name='variant_sample_information'),

    path('variant/<int:variant_id>/tag/<tag>/detail', views.variant_tag_detail, name='variant_tag_detail'),

    # Grids
    path('wiki/datatable/', DatabaseTableView.as_view(column_class=VariantWikiColumns),
         name='variant_wiki_datatable'),
    path('variant/<int:variant_id>/tag_counts/datatable/', DatabaseTableView.as_view(column_class=VariantTagCountsColumns),
         name='variant_tag_counts_datatable'),
    path('variant/<int:variant_id>/tag/<tag>/datatable/',
         DatabaseTableView.as_view(column_class=VariantTagDetailColumns),
         name='variant_tag_detail_datatable'),

    path('nearby/grid/<variant_id>/<genome_build_name>/<region_type>/<slug:op>/',
         JQGridView.as_view(grid=NearbyVariantsGrid, csv_download=True),
         name='nearby_variants_grid'),
    path('nearby/grid/<variant_id>/<genome_build_name>/<region_type>/<gene_symbol>/<slug:op>/',
         JQGridView.as_view(grid=NearbyVariantsGrid, csv_download=True),
         name='nearby_gene_variants_grid'),
    path('all_variants/grid/<genome_build_name>/<slug:op>/', JQGridView.as_view(grid=AllVariantsGrid, csv_download=True),
         name='all_variants_grid'),
    path('tags/grid/<genome_build_name>/<slug:op>/',
         JQGridView.as_view(grid=VariantTagsGrid, delete_row=True), name='variant_tags_grid'),
    path('tagged_variants/grid/<genome_build_name>/<slug:op>/',
         JQGridView.as_view(grid=TaggedVariantGrid, delete_row=True), name='tagged_variant_grid'),

    # Grid export
    path('tags/export/<genome_build_name>/', views.variant_tags_export, name='variant_tags_export'),
    path('tagged_variants/export/<genome_build_name>/', views.tagged_variant_export, name='tagged_variant_export'),

]
