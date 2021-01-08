from analysis.grids import AnalysesGrid, AnalysesVariantTagsGrid, NodeColumnSummaryGrid, \
    KaromappingAnalysesGrid, AnalysisTemplatesGrid, AnalysisNodeIssuesGrid, TaggedVariantGrid
from analysis.views import views, views_json, views_grid, \
    views_karyomapping, views_autocomplete
from library.django_utils.jqgrid_view import JQGridView
from variantgrid.perm_path import perm_path

urlpatterns = [
    perm_path('analyses/list/', views.analysis_list, name='analyses'),
    perm_path('analysis_templates/', views.analysis_templates, name='analysis_templates'),
    perm_path('<int:analysis_id>/', views.view_analysis, name='analysis'),
    perm_path('<int:analysis_id>/<int:active_node_id>/', views.view_analysis, name='analysis_node'),
    perm_path('clone_analysis/<int:analysis_id>/', views_json.clone_analysis, name='clone_analysis'),
    perm_path('create_analysis_from_template/<genome_build_name>', views.create_analysis_from_template, name='create_analysis_from_template'),
    perm_path('trio_wizard/<int:cohort_id>/<int:sample1_id>/<int:sample2_id>/<int:sample3_id>/', views.trio_wizard, name='trio_wizard'),

    # Templates
    perm_path('analysis_template/<pk>/save/', views_json.analysis_template_save, name='analysis_template_save'),
    perm_path('analysis_template/<pk>/list/', views.analysis_templates_list, name='analysis_templates_list'),

    # Node editor
    perm_path('node/view/<int:analysis_version>/<int:node_id>/<int:node_version>/<slug:extra_filters>/', views.node_view, name='node_view'),
    perm_path('node_update/<int:node_id>/', views_json.NodeUpdate.as_view(), name='node_update'),
    perm_path('node_debug/<int:analysis_version>/<int:node_id>/<int:node_version>/<slug:extra_filters>/', views.node_debug, name='node_debug'),
    perm_path('node_doc/<int:node_id>/', views.node_doc, name='node_doc'),
    perm_path('node_load/<int:node_id>/', views.node_load, name='node_load'),
    perm_path('node_cancel_load/<int:node_id>/', views.node_cancel_load, name='node_cancel_load'),
    perm_path('node/column_summary/<int:analysis_version>/<int:node_id>/<int:node_version>/<slug:extra_filters>/<slug:grid_column_name>/<int:significant_figures>/', views.node_column_summary, name='node_column_summary'),
    perm_path('node/node_snp_matrix/<int:node_id>/<int:node_version>/<slug:conversion>/<int:significant_figures>/', views.node_snp_matrix, name='node_snp_matrix'),
    perm_path('node/graph/<int:analysis_version>/<int:node_id>/<int:node_version>/<slug:graph_type_id>/<slug:cmap>/', views.node_data_graph, name='node_data_graph'),
    perm_path('node/cohort_zygosity_filters/<int:cohort_node_id>/<int:cohort_id>/', views.cohort_zygosity_filters, name='cohort_zygosity_filters'),
    perm_path('node/vcf_locus_filters/<int:node_id>/<int:vcf_id>/', views.vcf_locus_filters, name='vcf_locus_filters'),
    perm_path('node/sample_vcf_locus_filters/<int:node_id>/<int:sample_id>/', views.sample_vcf_locus_filters, name='sample_vcf_locus_filters'),
    perm_path('node/cohort_vcf_locus_filters/<int:node_id>/<int:cohort_id>/', views.cohort_vcf_locus_filters, name='cohort_vcf_locus_filters'),
    perm_path('node/pedigree_vcf_locus_filters/<int:node_id>/<int:pedigree_id>/', views.pedigree_vcf_locus_filters, name='pedigree_vcf_locus_filters'),

    # Node JSON
    perm_path('node/<int:node_id>/data/', views_json.node_data, name='node_data'),
    perm_path('analysis/<int:analysis_id>/node/create/<node_type>/', views_json.node_create, name='node_create'),
    perm_path('analysis/<int:analysis_id>/nodes/copy/', views_json.nodes_copy, name='nodes_copy'),
    perm_path('analysis/<int:analysis_id>/nodes/delete/', views_json.nodes_delete, name='nodes_delete'),
    perm_path('analysis/<int:analysis_id>/nodes/status/', views_json.nodes_status, name='nodes_status'),

    perm_path('create_filter_child/<int:node_id>/', views_json.create_filter_child, name='create_filter_child'),
    perm_path('create_extra_filter_child/<int:node_id>/<slug:extra_filters>/', views_json.create_extra_filter_child, name='create_extra_filter_child'),
    perm_path('create_selected_child/<int:node_id>/', views_json.create_selected_child, name='create_selected_child'),

    perm_path('analysis/<int:analysis_id>/node_versions/', views_json.analysis_node_versions, name='analysis_node_versions'),
    perm_path('analysis/<int:analysis_id>/edit_and_grid/', views.analysis_editor_and_grid, name='analysis_editor_and_grid'),
    perm_path('analysis/<int:analysis_id>/edit_and_grid/stand_alone/', views.stand_alone_analysis_editor_and_grid, name='standalone_analysis_editor_and_grid'),

    perm_path('analysis/<int:analysis_id>/set_panel_size/', views_json.analysis_set_panel_size, name='analysis_set_panel_size'),
    perm_path('analysis/<int:analysis_id>/settings/lock', views.analysis_settings_lock, name='analysis_settings_lock'),
    perm_path('analysis/<int:analysis_id>/settings/', views.view_analysis_settings, name='analysis_settings'),
    perm_path('analysis/<int:analysis_id>/settings_details_tab/', views.analysis_settings_details_tab, name='analysis_settings_details_tab'),
    perm_path('analysis/<int:analysis_id>/settings_node_counts_tab/', views.analysis_settings_node_counts_tab, name='analysis_settings_node_counts_tab'),
    perm_path('analysis/<int:analysis_id>/settings_template_run_tab/', views.analysis_settings_template_run_tab,
              name='analysis_settings_template_run_tab'),
    perm_path('analysis/<int:analysis_id>/reload/', views_json.analysis_reload, name='analysis_reload'),
    perm_path('analysis/<int:analysis_id>/input_samples/', views.analysis_input_samples, name='analysis_input_samples'),

    perm_path('node_graph/<int:node_id>/<int:graph_type_id>/<slug:cmap>/', views.node_graph, name='node_graph'),
    perm_path('column_summary_boxplot/<int:node_id>/<label>/<slug:variant_column>/', views.column_summary_boxplot, name='column_summary_boxplot'),
    perm_path('analysis/<int:analysis_id>/set_variant_tag/', views_json.set_variant_tag, name='set_variant_tag'),
    perm_path('set_variant_selected/<int:node_id>/', views_json.set_variant_selected, name='set_variant_selected'),
    perm_path('create_classification_from_variant_tag/<int:analysis_id>/<int:sample_id>/<int:variant_tag_id>/<transcript_id>/', views_json.create_classification_from_variant_tag, name='create_classification_with_transcript_from_variant_tag'),
    perm_path('create_classification_from_variant_tag/<int:analysis_id>/<int:sample_id>/<int:variant_tag_id>/', views_json.create_classification_from_variant_tag, name='create_classification_from_variant_tag'),

    # Node Data (bottom right window)
    perm_path('node_data_grid/cfg/<int:analysis_version>/<int:node_id>/<int:node_version>/<slug:extra_filters>/', views.node_data_grid, name='node_data_grid'),
    perm_path('node_async_wait/<int:analysis_version>/<int:node_id>/<int:node_version>/<slug:extra_filters>/', views.node_async_wait, name='node_async_wait'),
    perm_path('node_errors/<int:analysis_version>/<int:node_id>/<int:node_version>/<slug:extra_filters>/', views.node_errors, name='node_errors'),
    perm_path('node_method_description/<int:node_id>/<int:node_version>', views.node_method_description, name='node_method_description'),

    # Analysis templates
    perm_path('templates/variable/<int:node_id>/', views_json.analysis_template_variable, name='analysis_template_variable'),

    # Grids
    perm_path('node_grid/export/', views_grid.node_grid_export, name='node_grid_export'),
    perm_path('node_grid/cfg/<int:analysis_version>/<int:node_id>/<int:node_version>/<slug:extra_filters>/', views_grid.NodeGridConfig.as_view(), name='node_grid_config'),
    perm_path('node_grid/handler/', views_grid.NodeGridHandler.as_view(), name='node_grid_handler'),

    perm_path('node_column_summary/grid/<int:node_id>/<int:node_version>/<slug:extra_filters>/<slug:variant_column>/<int:significant_figures>/<slug:op>/', JQGridView.as_view(grid=NodeColumnSummaryGrid, csv_download=True), name='node_column_summary_grid'),
    perm_path('analyses/grid/<slug:op>/', JQGridView.as_view(grid=AnalysesGrid, delete_row=True), name='analyses_grid'),
    perm_path('analysis_templates/grid/<slug:op>/', JQGridView.as_view(grid=AnalysisTemplatesGrid, delete_row=True), name='analysis_templates_grid'),
    perm_path('tags/grid/<genome_build_name>/<slug:op>/', JQGridView.as_view(grid=AnalysesVariantTagsGrid, delete_row=True), name='analyses_variant_tags_grid'),
    perm_path('tagged_variants/grid/<genome_build_name>/<slug:op>/',
              JQGridView.as_view(grid=TaggedVariantGrid, delete_row=True), name='analysis_tagged_variant_grid'),

    perm_path('analysis_issues/grid/<slug:op>/',
              JQGridView.as_view(grid=AnalysisNodeIssuesGrid), name='analysis_node_issues_grid'),
    perm_path('analysis_issues', views.view_analysis_issues, name='analysis_issues'),

    # Mutational Signature
    perm_path('view_mutational_signature/<int:pk>/', views.view_mutational_signature, name='view_mutational_signature'),

    # Tagging
    perm_path('analyses_variant_tags/', views.analyses_variant_tags, name='analyses_variant_tags'),
    perm_path('analyses_variant_tags/<genome_build_name>', views.analyses_variant_tags, name='genome_build_analyses_variant_tags'),

    # karyomapping
    perm_path('karyomapping/analyses/', views_karyomapping.karyomapping_analyses, name='karyomapping_analyses'),
    perm_path('karyomapping/create_for_trio/<int:trio_id>/', views_karyomapping.create_karyomapping_analysis_for_trio_id, name='create_karyomapping_analysis_for_trio'),
    perm_path('karyomapping/view_karyomapping_analysis/<int:pk>/', views_karyomapping.view_karyomapping_analysis, name='view_karyomapping_analysis'),
    perm_path('karyomapping/view_karyomapping_gene/<int:pk>/', views_karyomapping.view_karyomapping_gene, name='view_karyomapping_gene'),
    perm_path('karyomapping/download_karyomapping_gene_csv/<int:pk>/', views_karyomapping.download_karyomapping_gene_csv, name='download_karyomapping_gene_csv'),

    perm_path('karyomapping/analyses/grid/<slug:op>/', JQGridView.as_view(grid=KaromappingAnalysesGrid, delete_row=True), name='karyomapping_analyses_grid'),

    # Autocompletes
    perm_path('autocomplete/Analysis/', views_autocomplete.AnalysisAutocompleteView.as_view(), name='analysis_autocomplete'),
    perm_path('autocomplete/AnalysisTemplate/', views_autocomplete.AnalysisTemplateAutocompleteView.as_view(),
              name='analysis_template_autocomplete'),
]
