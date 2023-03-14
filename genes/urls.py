from rest_framework.urlpatterns import format_suffix_patterns

from genes import views, views_autocomplete, views_rest
from genes.grids import GeneListGenesGrid, GenesGrid, QCGeneCoverageGrid, \
    UncoveredGenesGrid, GeneSymbolVariantsGrid, GeneSymbolWikiColumns, \
    GeneListColumns, CanonicalTranscriptCollectionColumns, CanonicalTranscriptColumns
from library.django_utils.jqgrid_view import JQGridView
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('genes', views.genes, name='genes'),
    perm_path('genes/<genome_build_name>', views.genes, name='genome_build_genes'),
    perm_path('view_gene/<gene_id>', views.view_gene, name='view_gene'),
    perm_path('view_gene_symbol/<gene_symbol>', views.view_gene_symbol, name='view_gene_symbol'),
    perm_path('view_gene_symbol/<gene_symbol>/<genome_build_name>/classification', views.view_classifications, name='view_gene_symbol_classifications'),
    perm_path('view_gene_symbol/<gene_symbol>/<genome_build_name>/classifications_download', views.export_classifications_gene_symbol, name='view_gene_symbol_classifications_download'),
    perm_path('view_gene_symbol/<gene_symbol>/<genome_build_name>', views.view_gene_symbol, name='view_gene_symbol_genome_build'),
    perm_path('view_transcript/<transcript_id>', views.view_transcript, name='view_transcript'),
    perm_path('view_transcript_version/<transcript_id>/<int:version>', views.view_transcript_version, name='view_transcript_version'),
    perm_path('view_transcript_accession/<transcript_accession>', views.view_transcript_accession, name='view_transcript_accession'),
    perm_path('gene_symbol_info_tab/<gene_symbol>/tool_tips/<tool_tips>', views.gene_symbol_info_tab, name='gene_symbol_info_tab'),
    perm_path('gene_lists', views.gene_lists, name='gene_lists'),
    perm_path('gene_lists/create_custom', views.create_custom_gene_list, name='create_custom_gene_list'),
    perm_path('view_gene_list/<int:gene_list_id>', views.view_gene_list, name='view_gene_list'),
    perm_path('view_canonical_transcript_collection/<pk>', views.view_canonical_transcript_collection, name='view_canonical_transcript_collection'),
    perm_path('qc_coverage', views.qc_coverage, name='qc_coverage'),
    perm_path('qc_coverage/<genome_build_name>', views.qc_coverage, name='genome_build_qc_coverage'),
    perm_path('gene_coverage_collection_graphs/<genome_build_name>/<slug:gene_symbol>', views.gene_coverage_collection_graphs, name='gene_symbol_coverage_collection_graphs'),
    perm_path('qc_gene_list_coverage_graphs/<genome_build_name>/<int:gene_list_id>', views.qc_gene_list_coverage_graphs, name='qc_gene_list_coverage_graphs'),
    perm_path('gene_grid/<path:columns_from_url>', views.gene_grid, name='passed_gene_grid'),
    perm_path('gene_grid', views.gene_grid, name='gene_grid'),
    perm_path('canonical_transcripts', views.canonical_transcripts, name='canonical_transcripts'),
    perm_path('sample_gene_lists_tab/<int:sample_id>', views.sample_gene_lists_tab, name='sample_gene_lists_tab'),
    perm_path('wiki', views.gene_wiki, name='gene_wiki'),

    perm_path('hotspot_graph/gene_symbol/<genome_build_name>/<gene_symbol>',
              views.HotspotGraphView.as_view(), name='gene_symbol_hotspot_graph'),
    perm_path('hotspot_graph/gene_symbol/<genome_build_name>/<gene_symbol>/transcript/<transcript_accession>',
              views.HotspotGraphView.as_view(), name='gene_symbol_transcript_version_hotspot_graph'),
    perm_path('hotspot_graph/gene/<genome_build_name>/<gene_id>',
              views.HotspotGraphView.as_view(), name='gene_hotspot_graph'),
    perm_path('hotspot_graph/transcript/<genome_build_name>/<transcript_id>',
              views.HotspotGraphView.as_view(), name='transcript_hotspot_graph'),
    perm_path('hotspot_graph/classifications/gene_symbol/<genome_build_name>/<gene_symbol>',
              views.ClassificationsHotspotGraphView.as_view(), name='classifications_gene_symbol_hotspot_graph'),
    perm_path('hotspot_graph/classifications/gene_symbol/<genome_build_name>/<gene_symbol>/transcript/<transcript_accession>',
              views.ClassificationsHotspotGraphView.as_view(), name='classifications_gene_symbol_transcript_version_hotspot_graph'),
    perm_path('hotspot_graph/classifications/gene/<genome_build_name>/<gene_id>',
              views.ClassificationsHotspotGraphView.as_view(), name='classifications_gene_hotspot_graph'),
    perm_path('hotspot_graph/classifications/transcript/<genome_build_name>/<transcript_id>',
              views.ClassificationsHotspotGraphView.as_view(), name='classifications_transcript_hotspot_graph'),
    perm_path('hotspot_graph/cohort/<int:cohort_id>/<transcript_id>',
              views.CohortHotspotGraphView.as_view(), name='cohort_hotspot_graph'),
    perm_path('hotspot_graph/public', views.PublicRUNX1HotspotGraphView.as_view(), name='public_hotspot_graph'),

    # Grids
    perm_path('wiki/datatable', DatabaseTableView.as_view(column_class=GeneSymbolWikiColumns),
              name='gene_wiki_datatable'),
    perm_path('gene/grid/<gene_symbol>/<genome_build_name>/<slug:op>/', JQGridView.as_view(grid=GeneSymbolVariantsGrid), name='gene_symbol_variants_grid'),

    perm_path('gene_lists/datatable/', DatabaseTableView.as_view(column_class=GeneListColumns),
              name='gene_lists_datatable'),
    perm_path('gene_list_genes/grid/<int:gene_list_id>/<slug:op>/', JQGridView.as_view(grid=GeneListGenesGrid, delete_row=True, csv_download=True), name='gene_list_genes_grid'),


    perm_path('canonical_transcript_collections/datatable/',
              DatabaseTableView.as_view(column_class=CanonicalTranscriptCollectionColumns),
              name='canonical_transcript_collections_datatable'),
    perm_path('canonical_transcript_collection/datatable/',
              DatabaseTableView.as_view(column_class=CanonicalTranscriptColumns),
              name='canonical_transcript_datatable'),
    perm_path('genes/grid/<genome_build_name>/<slug:op>/', JQGridView.as_view(grid=GenesGrid, csv_download=True), name='genes_grid'),
    perm_path('gene_coverage/grid/<int:gene_coverage_collection_id>/<slug:op>/', JQGridView.as_view(grid=QCGeneCoverageGrid), name='gene_coverage_collection_grid'),
    perm_path('gene_coverage/grid/<int:gene_coverage_collection_id>/<slug:op>/<path:gene_list_id_list>/',
              JQGridView.as_view(grid=QCGeneCoverageGrid), name='gene_coverage_collection_gene_list_grid'),
    perm_path('uncovered_genes/grid/<int:gene_coverage_collection_id>/<slug:op>/<path:gene_list_id_list>/min_depth/<int:min_depth>/', JQGridView.as_view(grid=UncoveredGenesGrid), name='uncovered_genes_grid'),

    perm_path('autocomplete/PanelAppPanel/aus', views_autocomplete.PanelAppPanelAusAutocompleteView.as_view(),
              name='panel_app_panel_aus_autocomplete'),
    perm_path('autocomplete/PanelAppPanel/eng', views_autocomplete.PanelAppPanelEngAutocompleteView.as_view(),
              name='panel_app_panel_eng_autocomplete'),
    perm_path('autocomplete/PanelAppPanel/', views_autocomplete.PanelAppPanelAutocompleteView.as_view(), name='panel_app_panel_autocomplete'),
    perm_path('autocomplete/Gene/', views_autocomplete.GeneAutocompleteView.as_view(), name='gene_autocomplete'),
    perm_path('autocomplete/Transcript/', views_autocomplete.TranscriptAutocompleteView.as_view(), name='transcript_autocomplete'),
    perm_path('autocomplete/GeneSymbol/', views_autocomplete.GeneSymbolAutocompleteView.as_view(), name='gene_symbol_autocomplete'),
    perm_path('autocomplete/GeneAnnotationRelease/', views_autocomplete.GeneAnnotationReleaseAutocompleteView.as_view(),
              name='gene_annotation_release_autocomplete'),
    perm_path('autocomplete/GeneList/category/', views_autocomplete.CategoryGeneListAutocompleteView.as_view(), name='category_gene_list_autocomplete'),
    perm_path('autocomplete/GeneList/', views_autocomplete.GeneListAutocompleteView.as_view(), name='gene_list_autocomplete'),
]

rest_urlpatterns = [
    perm_path('api/gene_list/modify/<int:pk>', views_rest.ModifyGeneListView.as_view(), name='api_modify_gene_list'),
    perm_path('api/gene_list/create', views_rest.CreateGeneListView.as_view(), name='api_create_gene_list'),
    perm_path('api/gene_list/<pk>', views_rest.GeneListView.as_view(), name='api_view_gene_list'),
    perm_path('api/panel_app/<pk>/gene_list', views_rest.PanelAppGeneListView.as_view(), name='api_view_panel_app_gene_list'),
    perm_path('api/named_gene_list/<category__name>/<name>', views_rest.NamedGeneListView.as_view(), name='api_named_gene_list'),
    perm_path('api/panel_app/gene_evidence/<int:server_id>/<gene_symbol>', views_rest.PanelAppGeneEvidenceView.as_view(), name='api_panel_app_gene_evidence'),
    perm_path('api/gene/batch_info', views_rest.BatchGeneInfoView.as_view(), name='api_batch_gene_info'),
    perm_path('api/gene_annotation_release/<int:release_id>/batch', views_rest.BatchGeneIdentifierForReleaseView.as_view(), name='api_batch_gene_identifiers_for_release'),
    perm_path('api/gene/info/<gene_symbol>', views_rest.GeneInfoView.as_view(), name='api_gene_info'),
    perm_path('api/text_to_gene_list', views_rest.TextToGeneListView.as_view(), name='api_text_to_gene_list'),
    perm_path('api/gene_annotation_release/<int:pk>', views_rest.GeneAnnotationReleaseView.as_view(), name='api_gene_annotation_release'),
    perm_path('api/sample_gene_list/<int:pk>', views_rest.SampleGeneListView.as_view(), name='api_sample_gene_list'),

]
urlpatterns += format_suffix_patterns(rest_urlpatterns)
