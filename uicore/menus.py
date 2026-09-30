"""
Menus as data (#2007): the top bar and every sub-menu are declared once, in MENUS, and a page's sub-menu and
highlighted items follow from its url name - templates never choose a menu.

A url name belongs to a menu as an item (listed in the sub-menu), as one of an item's `pages` (a detail page that
highlights that item, e.g. view_liftover_run under Liftover) or as one of the menu's own `pages` (in the menu but
under no item, e.g. view_allele under Variants). Per-deployment differences stay in URLS_NAME_REGISTER
(variantgrid/perm_path.py:get_visible_url_names); `condition` is only for the few items that depend on a setting.

Entry points: current_menu(url_name) and the menu_bar_main / menu_bar_sub tags in uicore/templatetags/ui_menus.py.
"""
from collections.abc import Callable
from dataclasses import dataclass
from typing import Optional

from django.conf import settings

from variantgrid.perm_path import get_visible_url_names


@dataclass(frozen=True)
class MenuItem:
    url_name: str
    title: Optional[str] = None  # defaults to the url name in title case
    admin_only: bool = False
    icon: Optional[str] = None
    href: Optional[str] = None  # link to this rather than reversing url_name
    external: bool = False
    method: str = 'get'  # 'post' renders a hidden form, for views that change state
    css_class: str = ''
    pages: tuple[str, ...] = ()  # url names of detail pages that highlight this item
    condition: Optional[Callable[[], bool]] = None  # the item isn't in this menu at all when False

    @property
    def is_available(self) -> bool:
        return self.condition is None or self.condition()

    @property
    def is_visible(self) -> bool:
        return self.is_available and (bool(self.href) or get_visible_url_names()[self.url_name])

    @property
    def display_title(self) -> str:
        return self.title or self.url_name.replace('_', ' ').title()

    def owns(self, url_name: str) -> bool:
        return url_name == self.url_name or url_name in self.pages


@dataclass(frozen=True)
class Menu:
    key: str
    title: str
    url_name: Optional[str]  # the top bar link; None keeps the menu out of the top bar (Settings)
    items: tuple[MenuItem, ...]
    pages: tuple[str, ...] = ()  # url names shown with this sub-menu without highlighting an item
    top_bar_condition: Optional[Callable[[], bool]] = None
    footer_template: Optional[str] = None

    @property
    def in_top_bar(self) -> bool:
        if not self.url_name or (self.top_bar_condition and not self.top_bar_condition()):
            return False
        return get_visible_url_names()[self.url_name]

    def visible_items(self) -> list[MenuItem]:
        return [item for item in self.items if item.is_visible]

    def owns(self, url_name: str) -> bool:
        """ Ignores URLS_NAME_REGISTER: a superuser reaching a hidden page still gets its menu """
        return url_name in self.pages or any(item.owns(url_name) for item in self.items if item.is_available)


def _seqauto_disabled() -> bool:
    return not settings.SEQAUTO_ENABLED


def _variants_menu_hidden() -> bool:
    """ Shariant has no Variants menu, so Liftover moves to Settings """
    return not get_visible_url_names()['variants']


LIFTOVER = MenuItem('liftover_runs', title="Liftover", admin_only=True, pages=('view_liftover_run',))
SEQUENCING_SOFTWARE_VERSIONS_PAGES = ('view_aligner', 'view_assay', 'view_library', 'view_sequencer',
                                      'view_variant_caller', 'view_variant_calling_pipeline')

MENUS: tuple[Menu, ...] = (
    Menu('sequencing', 'Sequencing', 'sequencing_data', items=(
        MenuItem('qc_coverage', title="Coverage", pages=('genome_build_qc_coverage',)),
        MenuItem('enrichment_kits_list', title="Enrichment Kits",
                 pages=('view_enrichment_kit', 'view_enrichment_kit_gene_coverage', 'view_gold_coverage_summary')),
        MenuItem('qc_data', title="QC Data", pages=('view_qc',)),
        MenuItem('qc_graphs', title="QC Graphs"),
        MenuItem('sequencing_data', title="Seq Runs",
                 pages=('view_sequencing_run', 'view_sequencing_run_tab', 'view_experiment', 'view_bam_file',
                        'view_unaligned_reads', 'view_single_sample_vcf', 'view_joint_called_vcf',
                        'view_tso500_pair')),
        MenuItem('sequencing_stats', title="Seq Stats", pages=('sequencing_stats_data',)),
        MenuItem('sequencing_software_versions', title="Seq / Software Versions",
                 pages=SEQUENCING_SOFTWARE_VERSIONS_PAGES),
    )),
    Menu('data', 'Data', 'data', items=(
        MenuItem('data'),
        MenuItem('upload', pages=('view_uploaded_file', 'view_upload_pipeline',
                                  'view_upload_pipeline_warnings_and_errors')),
    ), pages=('index', 'dashboard', 'search', 'staff_only', 'view_vcf', 'view_sample', 'bulk_group_permissions',
              'view_genome_build', 'view_contig', 'view_genomic_intervals', 'ontology_term',
              'mme_classification_panel', 'mme_view_submission', 'mme_view_inbound_match',
              'start_review', 'edit_review')),
    Menu('patients', 'Patients', 'patients', items=(
        MenuItem('cases', pages=('view_case',)),
        MenuItem('cohorts', pages=('view_cohort', 'cohort_gene_counts', 'cohort_hotspot')),
        MenuItem('patients', pages=('view_patient', 'patient_term_approvals', 'patient_term_approvals_offset',
                                    'patient_term_matches')),
        MenuItem('specimens', pages=('view_specimen',)),
        MenuItem('extractions', pages=('view_extraction', 'unmatched_extractions')),
        MenuItem('pedigrees', pages=('view_pedigree', 'view_ped_file')),
        MenuItem('trios', pages=('view_trio',)),
        MenuItem('quads', pages=('view_quad',)),
        MenuItem('duos', pages=('view_duo',)),
        MenuItem('patient_imports', pages=('view_patient_import', 'view_patient_record',
                                           'import_patient_records_details')),
    )),
    Menu('tests', 'Tests', 'pathology_tests', top_bar_condition=lambda: settings.PATHOLOGY_TESTS_ENABLED, items=(
        MenuItem('pathology_tests', title="Test Requests", pages=('view_pathology_test_order',)),
        MenuItem('manage_pathology_tests', title="Manage Tests",
                 pages=('view_pathology_test', 'view_pathology_test_version')),
    )),
    Menu('analysis', 'Analysis', 'analyses', items=(
        MenuItem('analysis_issues', title="Analysis Issues", admin_only=True),
        MenuItem('analyses'),
        MenuItem('analysis_templates', title="Templates",
                 pages=('analysis_templates_auto_launch', 'analysis_template_settings')),
        MenuItem('reanalysis_candidate_search', title="Candidate Search", pages=('new_reanalysis_candidate_search',)),
        MenuItem('karyomapping_analyses', title="Karyomapping", pages=('view_karyomapping_analysis',)),
    ), pages=('analysis', 'analysis_node', 'trio_wizard', 'quad_wizard', 'duo_wizard', 'view_mutational_signature')),
    Menu('classifications', 'Classifications', 'classifications', items=(
        MenuItem('activity', admin_only=True, pages=('activity_lab', 'activity_user', 'activity_discordance')),
        MenuItem('view_imported_allele_info', title="Imported Alleles", admin_only=True,
                 pages=('view_imported_allele_info_detail',)),
        MenuItem('clinvar_match', title="ClinVar Match", admin_only=True),
        MenuItem('classification_import_tool', title="Import Test", admin_only=True),
        MenuItem('hgvs_resolution_tool', title="HGVS Test", admin_only=True),
        MenuItem('condition_match_test', title="Condition Test", admin_only=True, pages=('condition_obsoletes',)),
        MenuItem('classification_view_metrics', title="View Metrics", admin_only=True),
        MenuItem('classification_reclassification_analytics', title="Reclassification", admin_only=True),
        MenuItem('classification_dashboard', title="Dashboard", icon="fas fa-home",
                 pages=('classification_dashboard_all',)),
        MenuItem('overlaps', title="Overlaps & Discordances",
                 pages=('overlap', 'overlap_history', 'action_overlap_review', 'discordance_report',
                        'discordance_report_deprecated', 'discordance_report_review_action')),
        MenuItem('classifications', pages=('view_classification', 'classification_history', 'classification_diff')),
        MenuItem('classification_candidate_search', title="Candidate Search",
                 pages=('new_classification_evidence_update_candidate_search',
                        'new_cross_sample_classification_candidate_search', 'view_candidate_search_run',
                        'classify_candidate')),
        MenuItem('clinvar_key_summary', title="ClinVar", pages=('clinvar_export',)),
        MenuItem('condition_matchings', title="Condition Matching",
                 pages=('condition_matching', 'condition_matchings_lab')),
        MenuItem('vus', title="VUS Resolution"),
        MenuItem('labs'),
        MenuItem('classification_graphs', title="Graphs / Stats"),
        MenuItem('classification_export', title="Export", pages=('classification_export_redirect',)),
        MenuItem('classification_grouping_export_config', title="Export NEW"),
        MenuItem('classification_upload_unmapped', title="Upload",
                 pages=('classification_upload_unmapped_lab', 'classification_upload_unmapped_status')),
    ), pages=('create_classification_for_variant', 'create_classification_for_variant_tag',
              'create_classification_from_hgvs', 'internal_lab_download', 'lab_gene_classification_counts')),
    Menu('genes', 'Genes', 'gene_lists', items=(
        MenuItem('genes', pages=('genome_build_genes',)),
        MenuItem('canonical_transcripts', pages=('view_canonical_transcript_collection',)),
        MenuItem('gene_grid', pages=('passed_gene_grid',)),
        MenuItem('gene_lists', pages=('view_gene_list',)),
        MenuItem('gene_wiki'),
    ), pages=('view_gene', 'view_gene_symbol', 'view_gene_symbol_genome_build', 'view_transcript',
              'view_transcript_version')),
    Menu('variants', 'Variants', 'variant_tags', items=(
        LIFTOVER,
        MenuItem('variant_tags', title="Tagged Variants", pages=('genome_build_variant_tags', 'tag_stats')),
        MenuItem('manual_variant_entry', title="Enter Variants", pages=('watch_manual_variant_entry',)),
        MenuItem('variants', pages=('genome_build_variants',)),
        MenuItem('variant_wiki', title="Variant Wiki", pages=('genome_build_variant_wiki',)),
    ), pages=('view_allele', 'view_allele_compact', 'view_variant', 'view_variant_genome_build',
              'view_variant_annotation_history', 'nearby_variants_annotation_version')),
    Menu('annotation', 'Annotation', 'annotation', items=(
        MenuItem('annotation'),
        MenuItem('variant_annotation_runs', title="Pipeline Runs", admin_only=True, pages=('view_annotation_run',)),
        MenuItem('view_annotation_descriptions', title="Descriptions",
                 pages=('view_annotation_descriptions_genome_build',)),
        MenuItem('pathogenicity_thresholds', title="Pathogenicity Thresholds"),
        MenuItem('annotation_versions', title="Versions"),
    ), pages=('view_citation',)),
    # Reached from the username in the navbar rather than the top bar
    Menu('settings', 'Settings', None, footer_template="uicore/menus/logout.html", items=(
        MenuItem('admin', href="/admin/", admin_only=True, external=True),
        MenuItem('eventlog', title="Event Log", admin_only=True),
        MenuItem('email_log', title="Email Log", admin_only=True, pages=('email_detail',)),
        MenuItem('keycloak_admin', title="Keycloak", admin_only=True),
        MenuItem('server_status', admin_only=True),
        MenuItem(LIFTOVER.url_name, title=LIFTOVER.title, admin_only=True, pages=LIFTOVER.pages,
                 condition=_variants_menu_hidden),
        MenuItem('change_password', pages=('password_change_done',)),
        MenuItem('custom_columns', pages=('view_custom_columns',)),
        MenuItem('tag_settings', pages=('view_tag_config_collection', 'tag_merge')),
        MenuItem('igv_integration', title="IGV Integration"),
        MenuItem('sequencing_software_versions', title="Sequencing / Software Versions",
                 pages=SEQUENCING_SOFTWARE_VERSIONS_PAGES, condition=_seqauto_disabled),
        MenuItem('view_user_settings', title="User Settings"),
        MenuItem('changelog'),
        MenuItem('version'),
    ), pages=('view_user', 'view_group', 'view_lab', 'view_organization', 'manual_migrations')),
    # Reached from the inbox icon in the navbar
    Menu('messages', 'Messages', None, items=(
        MenuItem('messages_inbox', title="Inbox", pages=('messages_detail',)),
        MenuItem('messages_trash', title="Trash"),
    ), pages=('messages_outbox', 'messages_compose', 'messages_compose_to')),
)


def current_menu(url_name: Optional[str]) -> Optional[Menu]:
    """ The menu that owns url_name, preferring one that is reachable (in the top bar, or has no top bar entry) -
        only Liftover and Seq / Software Versions are in two menus, and move to Settings when their own is off """
    if not url_name:
        return None
    owners = [menu for menu in MENUS if menu.owns(url_name)]
    for menu in owners:
        if menu.url_name is None or menu.in_top_bar:
            return menu
    return next(iter(owners), None)
