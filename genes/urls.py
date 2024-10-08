from rest_framework.urlpatterns import format_suffix_patterns

from genes.grids import GeneListGenesColumns, GenesGrid, QCGeneCoverageGrid, \
    UncoveredGenesGrid, GeneSymbolVariantsGrid, GeneSymbolWikiColumns, \
    GeneListColumns, CanonicalTranscriptCollectionColumns, CanonicalTranscriptColumns
from genes.views import views, views_autocomplete, views_rest
from genes.views.views_hotspot_graphs import HotspotGraphView, ClassificationsHotspotGraphView, CohortHotspotGraphView, \
    PublicRUNX1HotspotGraphView
from library.django_utils.jqgrid_view import JQGridView
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path

urlpatterns = [
    path('genes', views.genes, name='genes'),
    path('genes/<genome_build_name>', views.genes, name='genome_build_genes'),
    path('view_gene/<gene_id>', views.view_gene, name='view_gene'),
    path('view_gene_symbol/<gene_symbol>', views.view_gene_symbol, name='view_gene_symbol'),
    path('view_gene_symbol/<gene_symbol>/<genome_build_name>/classification', views.view_classifications, name='view_gene_symbol_classifications'),
    path('view_gene_symbol/<gene_symbol>/<genome_build_name>/classifications_download', views.export_classifications_gene_symbol, name='view_gene_symbol_classifications_download'),
    path('view_gene_symbol/<gene_symbol>/<genome_build_name>', views.view_gene_symbol, name='view_gene_symbol_genome_build'),
    path('view_transcript/<transcript_id>', views.view_transcript, name='view_transcript'),
    path('view_transcript_version/<transcript_id>/<int:version>', views.view_transcript_version, name='view_transcript_version'),
    path('view_transcript_accession/<transcript_accession>', views.view_transcript_accession, name='view_transcript_accession'),
    path('gene_symbol_info_tab/<gene_symbol>/tool_tips/<tool_tips>', views.gene_symbol_info_tab, name='gene_symbol_info_tab'),
    path('gene_lists', views.gene_lists, name='gene_lists'),
    path('gene_lists/create_custom', views.create_custom_gene_list, name='create_custom_gene_list'),
    path('view_gene_list/<int:gene_list_id>', views.view_gene_list, name='view_gene_list'),
    path('view_canonical_transcript_collection/<pk>', views.view_canonical_transcript_collection, name='view_canonical_transcript_collection'),
    path('qc_coverage', views.qc_coverage, name='qc_coverage'),
    path('qc_coverage/<genome_build_name>', views.qc_coverage, name='genome_build_qc_coverage'),
    path('gene_coverage_collection_graphs/<genome_build_name>/<slug:gene_symbol>', views.gene_coverage_collection_graphs, name='gene_symbol_coverage_collection_graphs'),
    path('qc_gene_list_coverage_graphs/<genome_build_name>/<int:gene_list_id>', views.qc_gene_list_coverage_graphs, name='qc_gene_list_coverage_graphs'),
    path('gene_grid/<path:columns_from_url>', views.gene_grid, name='passed_gene_grid'),
    path('gene_grid', views.gene_grid, name='gene_grid'),
    path('canonical_transcripts', views.canonical_transcripts, name='canonical_transcripts'),
    path('sample_gene_lists_tab/<int:sample_id>', views.sample_gene_lists_tab, name='sample_gene_lists_tab'),
    path('wiki', views.gene_wiki, name='gene_wiki'),

    path('hotspot_graph/gene_symbol/<genome_build_name>/<gene_symbol>',
         HotspotGraphView.as_view(), name='gene_symbol_hotspot_graph'),
    path('hotspot_graph/gene_symbol/<genome_build_name>/<gene_symbol>/transcript_accession/<transcript_accession>',
         HotspotGraphView.as_view(), name='gene_symbol_transcript_version_hotspot_graph'),
    path('hotspot_graph/gene/<genome_build_name>/<gene_id>',
         HotspotGraphView.as_view(), name='gene_hotspot_graph'),
    path('hotspot_graph/transcript/<genome_build_name>/<transcript_id>',
         HotspotGraphView.as_view(), name='transcript_hotspot_graph'),
    path('hotspot_graph/classifications/gene_symbol/<genome_build_name>/<gene_symbol>',
         ClassificationsHotspotGraphView.as_view(), name='classifications_gene_symbol_hotspot_graph'),
    path('hotspot_graph/classifications/gene_symbol/<genome_build_name>/<gene_symbol>/transcript_accession/<transcript_accession>',
         ClassificationsHotspotGraphView.as_view(), name='classifications_gene_symbol_transcript_version_hotspot_graph'),
    path('hotspot_graph/classifications/gene/<genome_build_name>/<gene_id>',
         ClassificationsHotspotGraphView.as_view(), name='classifications_gene_hotspot_graph'),
    path('hotspot_graph/classifications/transcript/<genome_build_name>/<transcript_id>',
         ClassificationsHotspotGraphView.as_view(), name='classifications_transcript_hotspot_graph'),
    path('hotspot_graph/cohort/<int:cohort_id>/<transcript_id>/<transcript_accession>',
         CohortHotspotGraphView.as_view(), name='cohort_transcript_version_hotspot_graph'),
    path('hotspot_graph/cohort/<int:cohort_id>/<transcript_id>',
         CohortHotspotGraphView.as_view(), name='cohort_transcript_hotspot_graph'),
    path('hotspot_graph/public', PublicRUNX1HotspotGraphView.as_view(), name='public_hotspot_graph'),

    # Grids
    path('wiki/datatable', DatabaseTableView.as_view(column_class=GeneSymbolWikiColumns),
         name='gene_wiki_datatable'),
    path('gene/grid/<gene_symbol>/<genome_build_name>/<slug:op>/', JQGridView.as_view(grid=GeneSymbolVariantsGrid), name='gene_symbol_variants_grid'),

    path('gene_lists/datatable/', DatabaseTableView.as_view(column_class=GeneListColumns),
         name='gene_lists_datatable'),
    path('gene_list_genes/grid/<int:gene_list_id>', DatabaseTableView.as_view(column_class=GeneListGenesColumns), name='gene_list_genes_datatable'),

    path('canonical_transcript_collections/datatable/',
         DatabaseTableView.as_view(column_class=CanonicalTranscriptCollectionColumns),
         name='canonical_transcript_collections_datatable'),
    path('canonical_transcript_collection/datatable/',
         DatabaseTableView.as_view(column_class=CanonicalTranscriptColumns),
         name='canonical_transcript_datatable'),
    path('genes/grid/<genome_build_name>/<slug:op>/', JQGridView.as_view(grid=GenesGrid, csv_download=True), name='genes_grid'),
    path('gene_coverage/grid/<int:gene_coverage_collection_id>/<slug:op>/', JQGridView.as_view(grid=QCGeneCoverageGrid), name='gene_coverage_collection_grid'),
    path('gene_coverage/grid/<int:gene_coverage_collection_id>/<slug:op>/<path:gene_list_id_list>/',
         JQGridView.as_view(grid=QCGeneCoverageGrid), name='gene_coverage_collection_gene_list_grid'),
    path('uncovered_genes/grid/<int:gene_coverage_collection_id>/<slug:op>/<path:gene_list_id_list>/min_depth/<int:min_depth>/', JQGridView.as_view(grid=UncoveredGenesGrid), name='uncovered_genes_grid'),

    path('autocomplete/PanelAppPanel/aus', views_autocomplete.PanelAppPanelAusAutocompleteView.as_view(),
         name='panel_app_panel_aus_autocomplete'),
    path('autocomplete/PanelAppPanel/eng', views_autocomplete.PanelAppPanelEngAutocompleteView.as_view(),
         name='panel_app_panel_eng_autocomplete'),
    path('autocomplete/PanelAppPanel/', views_autocomplete.PanelAppPanelAutocompleteView.as_view(), name='panel_app_panel_autocomplete'),
    path('autocomplete/Gene/', views_autocomplete.GeneAutocompleteView.as_view(), name='gene_autocomplete'),
    path('autocomplete/Transcript/', views_autocomplete.TranscriptAutocompleteView.as_view(), name='transcript_autocomplete'),
    path('autocomplete/GeneSymbol/', views_autocomplete.GeneSymbolAutocompleteView.as_view(), name='gene_symbol_autocomplete'),
    path('autocomplete/GeneAnnotationRelease/', views_autocomplete.GeneAnnotationReleaseAutocompleteView.as_view(),
         name='gene_annotation_release_autocomplete'),
    path('autocomplete/GeneList/category/', views_autocomplete.CategoryGeneListAutocompleteView.as_view(), name='category_gene_list_autocomplete'),
    path('autocomplete/GeneList/', views_autocomplete.GeneListAutocompleteView.as_view(), name='gene_list_autocomplete'),
]

rest_urlpatterns = [
    path('api/gene_list/modify/<int:pk>', views_rest.ModifyGeneListView.as_view(), name='api_modify_gene_list'),
    path('api/gene_list/create', views_rest.CreateGeneListView.as_view(), name='api_create_gene_list'),
    path('api/gene_list/<pk>', views_rest.GeneListView.as_view(), name='api_view_gene_list'),
    path('api/panel_app/<pk>/gene_list', views_rest.PanelAppGeneListView.as_view(), name='api_view_panel_app_gene_list'),
    path('api/named_gene_list/<category__name>/<name>', views_rest.NamedGeneListView.as_view(), name='api_named_gene_list'),
    path('api/panel_app/gene_evidence/<int:server_id>/<gene_symbol>', views_rest.PanelAppGeneEvidenceView.as_view(), name='api_panel_app_gene_evidence'),
    path('api/gene/batch_info', views_rest.BatchGeneInfoView.as_view(), name='api_batch_gene_info'),
    path('api/gene_annotation_release/<int:release_id>/batch', views_rest.BatchGeneIdentifierForReleaseView.as_view(), name='api_batch_gene_identifiers_for_release'),
    path('api/gene/info/<gene_symbol>', views_rest.GeneInfoView.as_view(), name='api_gene_info'),
    path('api/text_to_gene_list', views_rest.TextToGeneListView.as_view(), name='api_text_to_gene_list'),
    path('api/gene_annotation_release/<int:pk>', views_rest.GeneAnnotationReleaseView.as_view(), name='api_gene_annotation_release'),
    path('api/sample_gene_list/<int:pk>', views_rest.SampleGeneListView.as_view(), name='api_sample_gene_list'),

]
urlpatterns += format_suffix_patterns(rest_urlpatterns)
