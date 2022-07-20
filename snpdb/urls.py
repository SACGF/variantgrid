from django.contrib.auth.views import PasswordChangeView, PasswordChangeDoneView
from django.urls import include
from django.urls.conf import path
from rest_framework.documentation import include_docs_urls
from rest_framework.urlpatterns import format_suffix_patterns

from library.django_utils.jqgrid_view import JQGridView
from snpdb.grids import CohortListGrid, CohortSampleListGrid, SamplesListGrid, GenomicIntervalsListGrid, \
    CustomColumnsCollectionListGrid, TriosListGrid, VCFListGrid
from snpdb.views import views, views_json, views_rest, \
    views_autocomplete
from variantgrid.perm_path import perm_path

urlpatterns = [
    # public pages
    perm_path('public_global_sample_gene_matrix', views.public_global_sample_gene_matrix, name='public_global_sample_gene_matrix'),

    perm_path('data', views.data, name='data'),
    perm_path('vcfs', views.vcfs, name='vcfs'),
    perm_path('samples', views.samples, name='samples'),
    perm_path('bed_files', views.bed_files, name='bed_files'),
    perm_path('job_status/<job_id>', views_json.job_status, name='job_state'),
    perm_path('cached_generated_file/check/<cgf_id>', views_json.cached_generated_file_check, name='cached_generated_file_check'),
    perm_path('cached_generated_file/delete', views.cached_generated_file_delete, name='cached_generated_file_delete'),
    perm_path('cohort/create_cohort_genotype/<int:cohort_id>', views_json.create_cohort_genotype, name='create_cohort_genotype'),
    perm_path('cohort/create_sub_cohort/<int:cohort_id>', views_json.create_sub_cohort, name='create_sub_cohort'),
    perm_path('cohorts', views.cohorts, name='cohorts'),
    perm_path('view_vcf/<int:vcf_id>', views.view_vcf, name='view_vcf'),
    perm_path('vcf/<int:vcf_id>/populate_clingen_alleles', views_json.vcf_populate_clingen_alleles, name='vcf_populate_clingen_alleles'),
    perm_path('vcf/<int:vcf_id>/change_zygosity_count/<int:vzcc_id>/<operation>', views_json.vcf_change_zygosity_count, name='vcf_change_zygosity_count'),
    perm_path('get_patient_upload_csv_for_vcf/<int:pk>', views.get_patient_upload_csv_for_vcf, name='get_patient_upload_csv_for_vcf'),

    perm_path('view_sample/<int:sample_id>', views.view_sample, name='view_sample'),
    perm_path('view_genomic_intervals/<int:genomic_intervals_collection_id>', views.view_genomic_intervals, name='view_genomic_intervals'),
    perm_path('view_cohort_details_tab/<int:cohort_id>', views.view_cohort_details_tab, name='view_cohort_details_tab'),
    perm_path('view_cohort/<int:cohort_id>', views.view_cohort, name='view_cohort'),
    perm_path('cohort/hotspot/<int:cohort_id>', views.cohort_hotspot, name='cohort_hotspot'),
    perm_path('cohort/gene_counts/<int:cohort_id>', views.cohort_gene_counts, name='cohort_gene_counts'),
    perm_path('cohort/gene_counts/matrix/<int:cohort_id>/<gene_count_type_id>/<int:gene_list_id>', views.cohort_gene_counts_matrix, name='cohort_gene_counts_matrix'),
    perm_path('cohort/sort/<int:cohort_id>', views.cohort_sort, name='cohort_sort'),
    perm_path('cohort/sample_count/<int:cohort_id>', views_json.cohort_sample_count, name='cohort_sample_count'),
    perm_path('cohort_sample_edit/<int:cohort_id>', views.cohort_sample_edit, name='cohort_sample_edit'),
    perm_path('trios', views.trios, name='trios'),
    perm_path('view_trio/<int:pk>', views.view_trio, name='view_trio'),
    perm_path('sample_files_tab/<int:sample_id>', views.sample_files_tab, name='sample_files_tab'),
    perm_path('sample_variants_tab/<int:sample_id>', views.sample_variants_tab, name='sample_variants_tab'),
    perm_path('sample_variants_gene_detail/<int:sample_id>/<gene_symbol>', views.sample_variants_gene_detail, name='sample_variants_gene_detail'),
    perm_path('sample_graphs_tab/<int:sample_id>', views.sample_graphs_tab, name='sample_graphs_tab'),
    perm_path('sample_permissions_tab/<int:sample_id>', views.sample_permissions_tab, name='sample_permissions_tab'),
    perm_path('messages_bulk_delete', views.messages_bulk_delete, name='messages_bulk_delete'),
    perm_path('manual_variant_entry', views.manual_variant_entry, name='manual_variant_entry'),
    perm_path('watch_manual_variant_entry/<int:pk>', views.watch_manual_variant_entry, name='watch_manual_variant_entry'),
    perm_path('group_permissions/<class_name>/<int:primary_key>', views.group_permissions, name='group_permissions'),
    perm_path('group_permissions_object_delete/<class_name>/<int:primary_key>', views.group_permissions_object_delete, name='group_permissions_object_delete'),
    perm_path('group_permissions/bulk/<class_name>', views.bulk_group_permissions, name='bulk_group_permissions'),
    perm_path('settings/custom_columns/', views.custom_columns, name='custom_columns'),
    perm_path('settings/custom_columns/view/<int:custom_columns_collection_id>', views.view_custom_columns, name='view_custom_columns'),
    perm_path('settings/custom_columns/clone/<int:custom_columns_collection_id>', views_json.clone_custom_columns, name='clone_custom_columns'),
    perm_path('settings/set_user_row_config', views.set_user_row_config, name='set_user_row_config'),
    perm_path('settings/set_user_data_grid_config', views.set_user_data_grid_config, name='set_user_data_grid_config'),
    perm_path('settings/tags', views.tag_settings, name='tag_settings'),
    perm_path('settings/set_user_tag_color', views.set_user_tag_color, name='set_user_tag_color'),
    perm_path('settings/igv_integration', views.igv_integration, name='igv_integration'),
    perm_path('settings/password/password_change', PasswordChangeView.as_view(
        template_name='snpdb/settings/password_change.html',
    ), name='change_password'),
    perm_path('settings/password/password_change_done', PasswordChangeDoneView.as_view(template_name='snpdb/settings/password_change_done.html'), name='password_change_done'),
    perm_path('settings/user', views.view_user_settings, name='view_user_settings'),
    perm_path('settings/node_counts_tab/user', views.user_settings_node_counts_tab, name='user_settings_node_counts_tab'),
    perm_path('settings/node_counts_tab/lab/<int:pk>', views.lab_settings_node_counts_tab, name='lab_settings_node_counts_tab'),
    perm_path('settings/node_counts_tab/organization/<int:pk>', views.organization_settings_node_counts_tab, name='organization_settings_node_counts_tab'),

    perm_path('user/<pk>', views.view_user, name='view_user'),
    perm_path('lab/<int:lab_id>', views.view_lab, name='view_lab'),
    perm_path('clinvar_key/<str:pk>', views.view_clinvar_key, name='clinvar_key'),
    perm_path('organization/<int:organization_id>', views.view_organization, name='view_organization'),
    perm_path('help_static_page/<path:page_name>', views.help_static_page, name='help_static_page'),
    perm_path('genomic_intervals_graph/<int:genomic_intervals_collection_id>', views.genomic_intervals_graph, name='genomic_intervals_graph'),
    perm_path('chrom_density_graph/<int:sample_id>/<slug:cmap>', views.chrom_density_graph, name='chrom_density_graph'),
    perm_path('homozygosity_graph/<int:sample_id>/<slug:cmap>', views.homozygosity_graph, name='homozygosity_graph'),
    perm_path('sample_allele_frequency_histogram_graph/<int:sample_id>/<int:min_read_depth>', views.sample_allele_frequency_histogram_graph, name='sample_allele_frequency_histogram_graph'),
    perm_path('staff_only', views.staff_only, name='staff_only'),
    perm_path('tag_autocomplete_form', views.tag_autocomplete_form, name='tag_autocomplete_form'),
    perm_path('celery/wait_for_task/<celery_task>/<int:sleep_ms>/<path:redirect_url>', views.wait_for_task, name='wait_for_task'),
    perm_path('wiki_save/<class_name>/<unique_keyword>/<unique_value>', views.wiki_save, name='wiki_save'),
    perm_path('labs', views.labs, name='labs'),
    perm_path('labs_graph_detail', views.labs_graph_detail, name='labs_graph_detail'),
    perm_path('maps', views.maps, name='maps'),
    perm_path('index', views.index, name='index'),
    perm_path('user_global_sample_gene_matrix', views.user_global_sample_gene_matrix, name='user_global_sample_gene_matrix'),

    # Grids
    perm_path('cohort/grid/<slug:op>/', JQGridView.as_view(grid=CohortListGrid, delete_row=True), name='cohort_grid'),
    perm_path('cohort_sample/grid/<int:cohort_id>/<slug:op>/', JQGridView.as_view(grid=CohortSampleListGrid), name='cohort_sample_grid'),
    perm_path('sample/grid/<slug:op>/', JQGridView.as_view(grid=SamplesListGrid, delete_row=True), name='samples_grid'),
    perm_path('genomic_intervals/grid/<slug:op>/', JQGridView.as_view(grid=GenomicIntervalsListGrid, delete_row=True), name='genomic_intervals_grid'),
    perm_path('settings/custom_columns/grid/<slug:op>/', JQGridView.as_view(grid=CustomColumnsCollectionListGrid, delete_row=True), name='custom_columns_grid'),
    perm_path('trio/grid/<slug:op>/', JQGridView.as_view(grid=TriosListGrid, delete_row=True), name='trio_grid'),
    perm_path('vcfs/grid/<slug:op>/', JQGridView.as_view(grid=VCFListGrid, delete_row=True), name='vcfs_grid'),

    # Autocompletes
    perm_path('autocomplete/Cohort/', views_autocomplete.CohortAutocompleteView.as_view(), name='cohort_autocomplete'),
    perm_path('autocomplete/CustomColumn/', views_autocomplete.CustomColumnAutocompleteView.as_view(), name='custom_column_autocomplete'),
    perm_path('autocomplete/GenomicIntervalsCollection/', views_autocomplete.GenomicIntervalsCollectionAutocompleteView.as_view(), name='genomic_intervals_collection_autocomplete'),
    perm_path('autocomplete/Project/', views_autocomplete.ProjectAutocompleteView.as_view(), name='project_autocomplete'),
    perm_path('autocomplete/Sample/', views_autocomplete.SampleAutocompleteView.as_view(), name='sample_autocomplete'),
    perm_path('autocomplete/Tag/', views_autocomplete.TagAutocompleteView.as_view(), name='tag_autocomplete'),
    perm_path('autocomplete/Trio/', views_autocomplete.TrioAutocompleteView.as_view(), name='trio_autocomplete'),
    perm_path('autocomplete/User/', views_autocomplete.UserAutocompleteView.as_view(), name='user_autocomplete'),
    perm_path('autocomplete/Username/', views_autocomplete.UsernameAutocompleteView.as_view(), name='username_autocomplete'),
    perm_path('autocomplete/Lab/', views_autocomplete.LabAutocompleteView.as_view(), name='lab_autocomplete'),
    perm_path('autocomplete/VCF/', views_autocomplete.VCFAutocompleteView.as_view(), name='vcf_autocomplete'),

    # Debug dev help
    perm_path('ajax_hello_world/<str:data>', views.ajax_hello_world, name='ajax_hello_world'),
]

rest_urlpatterns = [
    path('api-auth/', include('rest_framework.urls')),
    perm_path('api/sample_variant_zygosity/<int:sample_id>/<int:variant_id>', views_rest.VariantZygosityForSampleView.as_view(), name='variant_zygosity_for_sample'),
    perm_path('api/trio/<pk>', views_rest.TrioView.as_view(), name='api_view_trio'),
    perm_path('api/variant_allele_for_variant/<int:variant_id>/<genome_build_name>',
              views_rest.VariantAlleleForVariantView.as_view(), name='variant_allele_for_variant'),
    perm_path('api/project/create', views_rest.ProjectViewSet.as_view({"post": "create"}), name='api_project_create'),
    perm_path('docs/', include_docs_urls(title='VariantGrid API', public=True, authentication_classes=[], permission_classes=[]), name='docs'),
]

urlpatterns += format_suffix_patterns(rest_urlpatterns)
