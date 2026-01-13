from django.contrib.auth.views import PasswordChangeView, PasswordChangeDoneView
from django.urls import include
from django.urls.conf import path as path_standard
from rest_framework.permissions import AllowAny
from rest_framework.renderers import OpenAPIRenderer
from rest_framework.schemas import get_schema_view

from library.django_utils.jqgrid_view import JQGridView
from library.preview_request import preview_view
from snpdb.grids import CohortListGrid, CohortSampleListGrid, SamplesListGrid, GenomicIntervalsListGrid, \
    CustomColumnsCollectionColumns, TriosListGrid, VCFListGrid, TagColorsCollectionColumns, \
    LiftoverRunColumns, LiftoverRunAlleleLiftoverColumns, AlleleLiftoverFailureColumns, \
    ManualVariantEntryCollectionColumns, SampleColumns
from snpdb.views import views, views_json, views_rest, views_autocomplete
from snpdb.views.datatable_view import DatabaseTableView
from variantgrid.perm_path import path

schema = get_schema_view(
    title="VariantGrid API",
    version="4.0.0",
    public=True,
    permission_classes=[AllowAny],
    renderer_classes=[OpenAPIRenderer],  # raw OpenAPI JSON
)


urlpatterns = [
    # public pages
    path('public_global_sample_gene_matrix', views.public_global_sample_gene_matrix, name='public_global_sample_gene_matrix'),

    path('data', views.data, name='data'),
    path('vcfs', views.vcfs, name='vcfs'),
    path('samples', views.samples, name='samples'),
    path('bed_files', views.bed_files, name='bed_files'),
    path('job_status/<job_id>', views_json.job_status, name='job_state'),
    path('cached_generated_file/check/<cgf_id>', views_json.cached_generated_file_check, name='cached_generated_file_check'),
    path('cached_generated_file/delete', views.cached_generated_file_delete, name='cached_generated_file_delete'),
    path('cohort/create_cohort_genotype/<int:cohort_id>', views_json.create_cohort_genotype, name='create_cohort_genotype'),
    path('cohort/create_sub_cohort/<int:cohort_id>', views_json.create_sub_cohort, name='create_sub_cohort'),
    path('cohorts', views.cohorts, name='cohorts'),
    path('view_vcf/<int:vcf_id>', views.view_vcf, name='view_vcf'),
    path('vcf/<int:vcf_id>/populate_clingen_alleles', views_json.vcf_populate_clingen_alleles, name='vcf_populate_clingen_alleles'),
    path('vcf/<int:vcf_id>/change_zygosity_count/<int:vzcc_id>/<operation>', views_json.vcf_change_zygosity_count, name='vcf_change_zygosity_count'),
    path('get_patient_upload_csv_for_vcf/<int:pk>', views.get_patient_upload_csv_for_vcf, name='get_patient_upload_csv_for_vcf'),

    path('view_sample/<int:sample_id>', views.view_sample, name='view_sample'),
    path('view_genomic_intervals/<int:genomic_intervals_collection_id>', views.view_genomic_intervals, name='view_genomic_intervals'),
    path('view_cohort_details_tab/<int:cohort_id>', views.view_cohort_details_tab, name='view_cohort_details_tab'),
    path('view_cohort/<int:cohort_id>', views.view_cohort, name='view_cohort'),
    path('cohort/hotspot/<int:cohort_id>', views.cohort_hotspot, name='cohort_hotspot'),
    path('cohort/gene_counts/<int:cohort_id>', views.cohort_gene_counts, name='cohort_gene_counts'),
    path('cohort/gene_counts/matrix/<int:cohort_id>/<gene_count_type_id>/<int:gene_list_id>', views.cohort_gene_counts_matrix, name='cohort_gene_counts_matrix'),
    path('cohort/sort/<int:cohort_id>', views.cohort_sort, name='cohort_sort'),
    path('cohort/sample_count/<int:cohort_id>', views_json.cohort_sample_count, name='cohort_sample_count'),
    path('cohort_sample_edit/<int:cohort_id>', views.cohort_sample_edit, name='cohort_sample_edit'),
    path('trios', views.trios, name='trios'),
    path('view_trio/<int:pk>', views.view_trio, name='view_trio'),
    path('sample_files_tab/<int:sample_id>', views.sample_files_tab, name='sample_files_tab'),
    path('sample_variants_tab/<int:sample_id>', views.sample_variants_tab, name='sample_variants_tab'),
    path('sample_variants_gene_detail/<int:sample_id>/<gene_symbol>', views.sample_variants_gene_detail, name='sample_variants_gene_detail'),
    path('sample_graphs_tab/<int:sample_id>', views.sample_graphs_tab, name='sample_graphs_tab'),
    path('sample_permissions_tab/<int:sample_id>', views.sample_permissions_tab, name='sample_permissions_tab'),
    path('messages_bulk_delete', views.messages_bulk_delete, name='messages_bulk_delete'),
    path('manual_variant_entry', views.manual_variant_entry, name='manual_variant_entry'),
    path('watch_manual_variant_entry/<int:pk>', views.watch_manual_variant_entry, name='watch_manual_variant_entry'),
    path('group_permissions/<class_name>/<int:primary_key>', views.group_permissions, name='group_permissions'),
    path('group_permissions_object_delete/<class_name>/<int:primary_key>', views.group_permissions_object_delete, name='group_permissions_object_delete'),
    path('group_permissions/bulk/<class_name>', views.bulk_group_permissions, name='bulk_group_permissions'),
    path('settings/custom_columns/', views.custom_columns, name='custom_columns'),
    path('settings/custom_columns/view/<int:custom_columns_collection_id>', views.view_custom_columns, name='view_custom_columns'),
    path('settings/custom_columns/clone/<int:custom_columns_collection_id>', views_json.clone_custom_columns, name='clone_custom_columns'),
    path('settings/set_user_row_config', views.set_user_row_config, name='set_user_row_config'),
    path('settings/set_user_data_grid_config', views.set_user_data_grid_config, name='set_user_data_grid_config'),
    path('settings/tags', views.tag_settings, name='tag_settings'),
    path('settings/tags/collection/datatable', DatabaseTableView.as_view(column_class=TagColorsCollectionColumns),
         name='tag_color_collections_datatable'),
    path('settings/tags/collection/view/<int:tag_colors_collection_id>', views.view_tag_colors_collection,
         name='view_tag_colors_collection'),
    path('settings/tags/collection/clone/<int:tag_colors_collection_id>', views_json.clone_tag_colors_collection,
         name='clone_tag_colors_collection'),
    path('settings/tags/collection/set_tag_color/<int:tag_colors_collection_id>', views.set_tag_color, name='set_tag_color'),
    path('settings/igv_integration', views.igv_integration, name='igv_integration'),
    path('settings/password/password_change', PasswordChangeView.as_view(
        template_name='snpdb/settings/password_change.html',
    ), name='change_password'),
    path('settings/password/password_change_done', PasswordChangeDoneView.as_view(template_name='snpdb/settings/password_change_done.html'), name='password_change_done'),
    path('settings/user', views.view_user_settings, name='view_user_settings'),
    path('settings/node_counts_tab/user', views.user_settings_node_counts_tab, name='user_settings_node_counts_tab'),
    path('settings/node_counts_tab/lab/<int:pk>', views.lab_settings_node_counts_tab, name='lab_settings_node_counts_tab'),
    path('settings/node_counts_tab/organization/<int:pk>', views.organization_settings_node_counts_tab, name='organization_settings_node_counts_tab'),

    path('user/<pk>', views.view_user, name='view_user'),
    path('group/<pk>', views.view_group, name='view_group'),
    path('lab/<int:lab_id>', views.view_lab, name='view_lab'),
    path('clinvar_key/<str:pk>', views.view_clinvar_key, name='clinvar_key'),
    path('organization/<int:organization_id>', views.view_organization, name='view_organization'),
    path('help_static_page/<path:page_name>', views.help_static_page, name='help_static_page'),
    path('genomic_intervals_graph/<int:genomic_intervals_collection_id>', views.genomic_intervals_graph, name='genomic_intervals_graph'),
    path('chrom_density_graph/<int:sample_id>/<slug:cmap>', views.chrom_density_graph, name='chrom_density_graph'),
    path('homozygosity_graph/<int:sample_id>/<slug:cmap>', views.homozygosity_graph, name='homozygosity_graph'),
    path('sample_allele_frequency_histogram_graph/<int:sample_id>/<int:min_read_depth>', views.sample_allele_frequency_histogram_graph, name='sample_allele_frequency_histogram_graph'),
    path('staff_only', views.staff_only, name='staff_only'),
    path('tag_autocomplete_form', views.tag_autocomplete_form, name='tag_autocomplete_form'),
    path('celery/wait_for_task/<celery_task>/<int:sleep_ms>/<path:redirect_url>', views.wait_for_task, name='wait_for_task'),
    path('wiki_save/<class_name>/<unique_keyword>/<unique_value>', views.wiki_save, name='wiki_save'),
    path('labs', views.labs, name='labs'),
    path('labs_graph_detail', views.labs_graph_detail, name='labs_graph_detail'),
    path('liftover', views.liftover_runs, name='liftover_runs'),
    path('liftover/view_liftover_run/<liftover_run_id>', views.view_liftover_run, name='view_liftover_run'),
    path('maps', views.maps, name='maps'),
    path('index', views.index, name='index'),
    path('user_global_sample_gene_matrix', views.user_global_sample_gene_matrix, name='user_global_sample_gene_matrix'),

    path('view_genome_build/<genome_build_name>', views.view_genome_build, name='view_genome_build'),
    path('view_contig/<contig_accession>', views.view_contig, name='view_contig'),

    path('manual_variant_entry_collections/detail/<int:pk>', views.manual_variant_entry_collection_detail,
         name='manual_variant_entry_collection_detail'),

    # Grids
    path('cohort/grid/<slug:op>/', JQGridView.as_view(grid=CohortListGrid, delete_row=True), name='cohort_grid'),
    path('cohort_sample/grid/<int:cohort_id>/<slug:op>/', JQGridView.as_view(grid=CohortSampleListGrid), name='cohort_sample_grid'),
    path('sample/grid/<slug:op>/', JQGridView.as_view(grid=SamplesListGrid, delete_row=True), name='samples_grid'),
    path('genomic_intervals/grid/<slug:op>/', JQGridView.as_view(grid=GenomicIntervalsListGrid, delete_row=True), name='genomic_intervals_grid'),
    path('liftover/liftover_runs/datatable', DatabaseTableView.as_view(column_class=LiftoverRunColumns),
         name='liftover_runs_datatable'),
    path('liftover/allele_liftover/datatable', DatabaseTableView.as_view(column_class=LiftoverRunAlleleLiftoverColumns),
         name='allele_liftover_datatable'),
    path('liftover/allele_liftover_failures/datatable',
         DatabaseTableView.as_view(column_class=AlleleLiftoverFailureColumns),
         name='allele_liftover_failures_datatable'),
    path('manual_variant_entry_collections/datatable',
         DatabaseTableView.as_view(column_class=ManualVariantEntryCollectionColumns),
         name='manual_variant_entry_collections_datatable'),
    path('settings/custom_columns/collection/datatable',
         DatabaseTableView.as_view(column_class=CustomColumnsCollectionColumns),
         name='custom_columns_collections_datatable'),
    path('samples/datatable/',
         DatabaseTableView.as_view(column_class=SampleColumns),
         name='samples_datatable'),
    path('trio/grid/<slug:op>/', JQGridView.as_view(grid=TriosListGrid, delete_row=True), name='trio_grid'),
    path('vcfs/grid/<slug:op>/', JQGridView.as_view(grid=VCFListGrid, delete_row=True), name='vcfs_grid'),

    # Autocompletes
    path('autocomplete/Cohort/', views_autocomplete.CohortAutocompleteView.as_view(), name='cohort_autocomplete'),
    path('autocomplete/CustomColumn/', views_autocomplete.CustomColumnAutocompleteView.as_view(), name='custom_column_autocomplete'),
    path('autocomplete/GenomicIntervalsCollection/', views_autocomplete.GenomicIntervalsCollectionAutocompleteView.as_view(), name='genomic_intervals_collection_autocomplete'),
    path('autocomplete/Project/', views_autocomplete.ProjectAutocompleteView.as_view(), name='project_autocomplete'),
    path('autocomplete/Sample/', views_autocomplete.SampleAutocompleteView.as_view(), name='sample_autocomplete'),
    path('autocomplete/Tag/', views_autocomplete.TagAutocompleteView.as_view(), name='tag_autocomplete'),
    path('autocomplete/Trio/', views_autocomplete.TrioAutocompleteView.as_view(), name='trio_autocomplete'),
    path('autocomplete/User/', views_autocomplete.UserAutocompleteView.as_view(), name='user_autocomplete'),
    path('autocomplete/Username/', views_autocomplete.UsernameAutocompleteView.as_view(), name='username_autocomplete'),
    path('autocomplete/Lab/', views_autocomplete.LabAutocompleteView.as_view(), name='lab_autocomplete'),
    path('autocomplete/VCF/', views_autocomplete.VCFAutocompleteView.as_view(), name='vcf_autocomplete'),

    # Previews
    path('preview/<str:db>/<str:idx>', preview_view, name='preview_data'),

    # For Uptime Robot
    path('uptime_check', views.view_uptime, name='uptime_check'),

    # Debug dev help
    path('ajax_hello_world/<str:data>', views.ajax_hello_world, name='ajax_hello_world'),

    path_standard('api-auth/', include('rest_framework.urls')),
    path('api/sample_variant_zygosity/<int:sample_id>/<int:variant_id>', views_rest.VariantZygosityForSampleView.as_view(), name='variant_zygosity_for_sample'),
    path('api/trio/<pk>', views_rest.TrioView.as_view(), name='api_view_trio'),
    path('api/variant_allele_for_variant/<int:variant_id>/<genome_build_name>',
         views_rest.VariantAlleleForVariantView.as_view(), name='variant_allele_for_variant'),
    path('api/project/create', views_rest.ProjectViewSet.as_view({"post": "create"}), name='api_project_create'),
    path("docs/", schema, name="openapi-schema"),
]
