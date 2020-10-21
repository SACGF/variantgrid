from library.django_utils.jqgrid_view import JQGridView
from variantgrid.perm_path import perm_path
from variantopedia import views
from variantopedia.grids import VariantWikiGrid, AllVariantsGrid, NearbyVariantsGrid

urlpatterns = [
    perm_path('variants', views.variants, name='variants'),
    perm_path('dashboard', views.dashboard, name='dashboard'),
    perm_path('server_status', views.server_status, name='server_status'),
    perm_path('database_statistics', views.database_statistics, name='database_statistics'),
    perm_path('tagged', views.tagged, name='variantopedia_tagged'),
    perm_path('search', views.search, name='search'),
    perm_path('view_variant/<int:variant_id>', views.view_variant, name='view_variant'),
    perm_path('view_variant_annotation_history/<int:variant_id>', views.view_variant_annotation_history,
              name='view_variant_annotation_history'),
    perm_path('wiki', views.wiki, name='variantopedia_wiki'),
    perm_path('view_allele_from_variant/<int:variant_id>', views.view_allele_from_variant,
              name='view_allele_from_variant'),
    perm_path('view_allele/<int:pk>', views.view_allele, name='view_allele'),
    perm_path('view/<int:variant_id>', views.variant_details, name='variant_details'),
    perm_path('view/<int:variant_id>/<int:annotation_version_id>', views.variant_details_annotation_version,
              name='variant_details_annotation_version'),
    perm_path('nearby/<int:variant_id>', views.nearby_variants, name='nearby_variants'),
    perm_path('view/<int:variant_id>/<int:annotation_version_id>', views.nearby_variants,
              name='variant_details_annotation_version'),
    perm_path('gene_coverage/<slug:gene_symbol_id>', views.gene_coverage, name='gene_coverage'),
    perm_path('variant_sample_information/<int:variant_id>', views.variant_sample_information,
              name='variant_sample_information'),

    # Grids
    perm_path('wiki/grid/<slug:op>/', JQGridView.as_view(grid=VariantWikiGrid, csv_download=True),
              name='variantopedia_wiki_grid'),
    perm_path('nearby/grid/<variant_id>/<region_type>/<slug:op>/',
              JQGridView.as_view(grid=NearbyVariantsGrid,csv_download=True),
              name='nearby_variants_grid'),
    perm_path('all_variants/grid/<slug:op>/', JQGridView.as_view(grid=AllVariantsGrid, csv_download=True),
              name='all_variants_grid'),
]
