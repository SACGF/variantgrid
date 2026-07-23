from beacon import views, views_rest
from variantgrid.perm_path import path

urlpatterns = [
    # Framework (identity / config / map)
    path('', views_rest.BeaconInfoView.as_view(), name='beacon_index'),
    path('info', views_rest.BeaconInfoView.as_view(), name='beacon_info'),
    path('service-info', views_rest.BeaconServiceInfoView.as_view(), name='beacon_service_info'),
    path('configuration', views_rest.BeaconConfigurationView.as_view(), name='beacon_configuration'),
    path('entry_types', views_rest.BeaconEntryTypesView.as_view(), name='beacon_entry_types'),
    path('filtering_terms', views_rest.BeaconFilteringTermsView.as_view(), name='beacon_filtering_terms'),
    path('map', views_rest.BeaconMapView.as_view(), name='beacon_map'),

    # Model - genomic variants
    path('g_variants', views_rest.BeaconGVariantsView.as_view(), name='beacon_g_variants'),
    path('g_variants/<int:variant_id>', views_rest.BeaconGVariantByIdView.as_view(),
         name='beacon_g_variant_by_id'),

    # Outbound (§9.4): variant-page external-Beacon section
    path('external_beacons/<int:variant_id>/<genome_build_name>', views.external_beacons,
         name='external_beacons'),
]
