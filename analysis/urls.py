from analysis.grids import AnalysesGrid, NodeColumnSummaryGrid, KaromappingAnalysesGrid, AnalysisTemplatesGrid, \
    AnalysisNodeIssuesGrid, NodeOntologyGenesGrid, NodeGeneDiseaseClassificationGenesGrid, \
    NodeTissueExpressionGenesGrid, NodeTissueUniProtTissueSpecificityGenesGrid, NodeGeneListGenesColumns
from analysis.views import views, views_json, views_grid, views_karyomapping, views_autocomplete
from library.django_utils.jqgrid_view import JQGridView
from snpdb.views.datatable_view import DatabaseTableView
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
    perm_path('analysis_template/<pk>/settings/', views.analysis_template_settings, name='analysis_template_settings'),
    perm_path('analysis_template/<pk>/clone/', views_json.analysis_template_clone, name='analysis_template_clone'),


    # For node views below - a Node contains the analysis ID - so we don't need to pass analysis_id, but do so
    # to make it easier to debug errors if nodes were deleted

    # Node editors
    perm_path('<int:analysis_id>/<int:analysis_version>/node/view/<int:node_id>/<int:node_version>/<slug:extra_filters>/', views.node_view, name='node_view'),
    perm_path('<int:analysis_id>/node_update/<int:node_id>/', views_json.NodeUpdate.as_view(), name='node_update'),
    perm_path('<int:analysis_id>/<int:analysis_version>/node_debug/<int:node_id>/<int:node_version>/<slug:extra_filters>/', views.node_debug, name='node_debug'),
    perm_path('<int:analysis_id>/node_doc/<int:node_id>/', views.node_doc, name='node_doc'),
    perm_path('<int:analysis_id>/node_load/<int:node_id>/', views.node_load, name='node_load'),
    perm_path('<int:analysis_id>/node_cancel_load/<int:node_id>/', views.node_cancel_load, name='node_cancel_load'),
    perm_path('<int:analysis_id>/<int:analysis_version>/node/column_summary/<int:node_id>/<int:node_version>/<slug:extra_filters>/<slug:grid_column_name>/<int:significant_figures>/', views.node_column_summary, name='node_column_summary'),
    perm_path('<int:analysis_id>/node/node_snp_matrix/<int:node_id>/<int:node_version>/<slug:conversion>/<int:significant_figures>/', views.node_snp_matrix, name='node_snp_matrix'),
    perm_path('<int:analysis_id>/<int:analysis_version>/node/graph/<int:node_id>/<int:node_version>/<slug:graph_type_id>/<slug:cmap>/', views.node_data_graph, name='node_data_graph'),
    perm_path('<int:analysis_id>/node/cohort_zygosity_filters/<int:node_id>/<int:cohort_id>/', views.cohort_zygosity_filters, name='cohort_zygosity_filters'),
    perm_path('<int:analysis_id>/node/vcf_locus_filters/<int:node_id>/<int:vcf_id>/', views.vcf_locus_filters, name='vcf_locus_filters'),
    perm_path('<int:analysis_id>/node/sample_vcf_locus_filters/<int:node_id>/<int:sample_id>/', views.sample_vcf_locus_filters, name='sample_vcf_locus_filters'),
    perm_path('<int:analysis_id>/node/cohort_vcf_locus_filters/<int:node_id>/<int:cohort_id>/', views.cohort_vcf_locus_filters, name='cohort_vcf_locus_filters'),
    perm_path('<int:analysis_id>/node/pedigree_vcf_locus_filters/<int:node_id>/<int:pedigree_id>/', views.pedigree_vcf_locus_filters, name='pedigree_vcf_locus_filters'),

    # Node JSON
    perm_path('<int:analysis_id>/node/<int:node_id>/data/', views_json.node_data, name='node_data'),
    perm_path('<int:analysis_id>/node/create/<node_type>/', views_json.node_create, name='node_create'),
    perm_path('<int:analysis_id>/nodes/copy/', views_json.nodes_copy, name='nodes_copy'),
    perm_path('<int:analysis_id>/nodes/delete/', views_json.nodes_delete, name='nodes_delete'),
    perm_path('<int:analysis_id>/nodes/status/', views_json.nodes_status, name='nodes_status'),
    perm_path('<int:analysis_id>/node/tasks/', views_json.nodes_tasks, name='nodes_tasks'),

    perm_path('<int:analysis_id>/create_filter_child/<int:node_id>/', views_json.create_filter_child, name='create_filter_child'),
    perm_path('<int:analysis_id>/create_extra_filter_child/<int:node_id>/<slug:extra_filters>/', views_json.create_extra_filter_child, name='create_extra_filter_child'),
    perm_path('<int:analysis_id>/create_selected_child/<int:node_id>/', views_json.create_selected_child, name='create_selected_child'),

    perm_path('<int:analysis_id>/node_versions/', views_json.analysis_node_versions, name='analysis_node_versions'),
    perm_path('<int:analysis_id>/edit_and_grid/', views.analysis_editor_and_grid, name='analysis_editor_and_grid'),
    perm_path('<int:analysis_id>/edit_and_grid/stand_alone/', views.stand_alone_analysis_editor_and_grid, name='standalone_analysis_editor_and_grid'),

    perm_path('<int:analysis_id>/set_panel_size/', views_json.analysis_set_panel_size, name='analysis_set_panel_size'),
    perm_path('<int:analysis_id>/node_populate_clingen_alleles/<int:node_id>/', views_json.node_populate_clingen_alleles,
              name='node_populate_clingen_alleles'),
    perm_path('<int:analysis_id>/settings/lock', views.analysis_settings_lock, name='analysis_settings_lock'),
    perm_path('<int:analysis_id>/settings/', views.view_analysis_settings, name='analysis_settings'),
    perm_path('<int:analysis_id>/settings_details_tab/', views.analysis_settings_details_tab, name='analysis_settings_details_tab'),
    perm_path('<int:analysis_id>/settings_node_counts_tab/', views.analysis_settings_node_counts_tab, name='analysis_settings_node_counts_tab'),
    perm_path('<int:analysis_id>/settings_template_run_tab/', views.analysis_settings_template_run_tab,
              name='analysis_settings_template_run_tab'),
    perm_path('<int:analysis_id>/reload/', views_json.analysis_reload, name='analysis_reload'),
    perm_path('<int:analysis_id>/input_samples/', views.analysis_input_samples, name='analysis_input_samples'),
    perm_path('sample_patient_gene_disease/<int:sample_id>', views_json.sample_patient_gene_disease, name='sample_patient_gene_disease'),

    perm_path('<int:analysis_id>/node_graph/<int:node_id>/<int:graph_type_id>/<slug:cmap>/', views.node_graph, name='node_graph'),
    perm_path('<int:analysis_id>/column_summary_boxplot/<int:node_id>/<label>/<slug:variant_column>/', views.column_summary_boxplot, name='column_summary_boxplot'),
    perm_path('<int:analysis_id>/set_variant_selected/<int:node_id>/', views_json.set_variant_selected, name='set_variant_selected'),

    perm_path('set_variant_tag/<slug:location>/', views_json.set_variant_tag, name='set_variant_tag'),

    perm_path('<int:analysis_id>/classification/create_for_variant_tag/<int:variant_tag_id>', views.CreateClassificationForVariantTagView.as_view(),
              name='create_classification_for_variant_tag'),
    perm_path('<int:analysis_id>/create_classification/',
              views.create_classification_for_analysis, name='create_classification_for_analysis'),

    # Node Data (bottom right window)
    perm_path('<int:analysis_id>/<int:analysis_version>/node_data_grid/cfg/<int:node_id>/<int:node_version>/<slug:extra_filters>/', views.node_data_grid, name='node_data_grid'),
    perm_path('<int:analysis_id>/<int:analysis_version>/node_async_wait/<int:node_id>/<int:node_version>/<slug:extra_filters>/', views.node_async_wait, name='node_async_wait'),
    perm_path('<int:analysis_id>/<int:analysis_version>/node_errors/<int:node_id>/<int:node_version>/<slug:extra_filters>/', views.node_errors, name='node_errors'),
    perm_path('<int:analysis_id>/node_method_description/<int:node_id>/<int:node_version>', views.node_method_description, name='node_method_description'),

    # Analysis templates
    perm_path('<int:analysis_id>/templates/variable/<int:node_id>/', views_json.analysis_template_variable, name='analysis_template_variable'),

    # Grids
    perm_path('<int:analysis_id>/node_grid/export/', views_grid.node_grid_export, name='node_grid_export'),
    perm_path('<int:analysis_id>/<int:analysis_version>/node_grid/cfg/<int:node_id>/<int:node_version>/<slug:extra_filters>/', views_grid.NodeGridConfig.as_view(), name='node_grid_config'),
    perm_path('<int:analysis_id>/node_grid/handler/', views_grid.NodeGridHandler.as_view(), name='node_grid_handler'),

    perm_path('<int:analysis_id>/node_column_summary/grid/<int:node_id>/<int:node_version>/<slug:extra_filters>/<slug:variant_column>/<int:significant_figures>/<slug:op>/', JQGridView.as_view(grid=NodeColumnSummaryGrid, csv_download=True), name='node_column_summary_grid'),
    perm_path('analyses/grid/<slug:op>/', JQGridView.as_view(grid=AnalysesGrid, delete_row=True), name='analyses_grid'),

    perm_path('analysis_templates/grid/<slug:op>/', JQGridView.as_view(grid=AnalysisTemplatesGrid, delete_row=True), name='analysis_templates_grid'),
    perm_path('analysis_issues/grid/<slug:op>/',
              JQGridView.as_view(grid=AnalysisNodeIssuesGrid), name='analysis_node_issues_grid'),

    perm_path('<int:analysis_id>/node/ontology/genes/grid/<int:node_id>/<int:version>/<slug:op>/',
              JQGridView.as_view(grid=NodeOntologyGenesGrid), name='node_ontology_genes_grid'),
    perm_path('<int:analysis_id>/node/gene_disease_classification/grid/<int:node_id>/<int:version>/<slug:op>/',
              JQGridView.as_view(grid=NodeGeneDiseaseClassificationGenesGrid),
              name='node_gene_disease_classification_genes_grid'),

    perm_path('<int:analysis_id>/node/tissue/genes/grid/<int:node_id>/<int:version>/<slug:op>/',
              JQGridView.as_view(grid=NodeTissueExpressionGenesGrid), name='node_tissue_expression_genes_grid'),
    perm_path('<int:analysis_id>/node/tissue_uniprot/genes/grid/<int:node_id>/<int:version>/<slug:op>/',
              JQGridView.as_view(grid=NodeTissueUniProtTissueSpecificityGenesGrid), name='node_tissue_uniprot_genes_grid'),

    perm_path('<int:analysis_id>/node/<int:node_id>/<int:version>/gene_list_genes/<int:gene_list_id>',
              DatabaseTableView.as_view(column_class=NodeGeneListGenesColumns),
              name='analysis_node_gene_list_genes_datatable'),

    perm_path('analysis_issues', views.view_analysis_issues, name='analysis_issues'),

    # Mutational Signature
    perm_path('view_mutational_signature/<int:pk>/', views.view_mutational_signature, name='view_mutational_signature'),

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
