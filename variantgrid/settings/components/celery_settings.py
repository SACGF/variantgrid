from kombu import Exchange, Queue

from variantgrid.settings.components.secret_settings import get_secret

# Using RabbitMQ as Redis broker was re-executing long tasks (eg variant annotation)
# For instructions on setting up proper users:
# @see https://docs.celeryproject.org/en/latest/getting-started/brokers/rabbitmq.html#id4
# rabbitmqctl add_user variantgrid secretpassword
# rabbitmqctl add_vhost variantgrid
# rabbitmqctl set_user_tags variantgrid vg_tag
# rabbitmqctl set_permissions -p variantgrid variantgrid ".*" ".*" ".*"
#
# Then in /etc/variantgrid/settings_config.json
#    "CELERY": {
#        "broker_url": "amqp://variantgrid:secretpassword@localhost/variantgrid"
#    },
# Or if you just want the warning to go away: "broker_url": "amqp://guest:guest@localhost"


# Default is guest:guest default account see: https://www.rabbitmq.com/access-control.html
CELERY_BROKER_URL = get_secret("CELERY.broker_url")

CELERY_TASK_DEFAULT_QUEUE = 'db_workers'
CELERY_TASK_QUEUES = (
    Queue('analysis_workers', Exchange('analysis_workers'), routing_key='analysis_workers'),
    Queue('annotation_workers', Exchange('annotation_workers'), routing_key='annotation_workers'),
    # db_workers - only use the database, ie don't write to web server filesystem (so it could be moved to diff machine)
    Queue('db_workers', Exchange('db_workers'), routing_key='db_workers'),
    Queue('web_workers', Exchange('web_workers'), routing_key='web_workers'),
    Queue('scheduling_single_worker', Exchange('scheduling_single_worker'), routing_key='scheduling_single_worker'),
    Queue('variant_id_single_worker', Exchange('variant_id_single_worker'), routing_key='variant_id_single_worker'),
)

ANALYSIS_WORKERS = {'queue': 'analysis_workers', 'routing_key': 'analysis_workers'}
ANNOTATION_WORKERS = {'queue': 'annotation_workers', 'routing_key': 'annotation_workers'}
DB_WORKERS = {'queue': 'db_workers', 'routing_key': 'db_workers'}  # Default Queue
# Require reading from filesystem (upload directory or images & reports)
WEB_WORKERS = {'queue': 'web_workers', 'routing_key': 'web_workers'}

# This is used for inserting new Variants/Loci
# 1 worker to avoid race conditions so we end up with 1 and only 1 variant for chrom/post/ref/alt
VARIANT_ID_SINGLE_WORKER = {'queue': 'variant_id_single_worker', 'routing_key': 'variant_id_single_worker'}

# 1 worker, use this to schedule tasks and avoid race conditions
SCHEDULING_SINGLE_WORKER = {'queue': 'scheduling_single_worker', 'routing_key': 'scheduling_single_worker'}

CELERY_WORKER_NAMES = ['annotation_workers', 'db_workers', 'web_workers',
                       'scheduling_single_worker', 'variant_id_single_worker']
CELERY_ANALYSIS_WORKER_NAMES = ['analysis_workers']

CELERY_TASK_ROUTES = {
    # Analysis
    'analysis.tasks.auto_analysis_tasks.auto_run_analyses_for_sample': ANALYSIS_WORKERS,
    'analysis.tasks.auto_analysis_tasks.auto_run_analyses_for_vcf': ANALYSIS_WORKERS,
    'analysis.tasks.karyomapping_tasks.create_genome_karyomapping_for_trio': ANALYSIS_WORKERS,
    "analysis.tasks.node_update_tasks.update_node_task": ANALYSIS_WORKERS,
    "analysis.tasks.node_update_tasks.node_cache_task": ANALYSIS_WORKERS,
    # Venn equivalent of node_cache_task - chained ahead of the Venn node's update. Keep it on the
    # analysis pool (not the default db_workers) so cache building scales with the rest of the
    # pipeline rather than starving on db_workers and triggering lease-expiry re-dispatch churn.
    "analysis.models.nodes.filters.venn_node.venn_cache_count": ANALYSIS_WORKERS,
    "analysis.tasks.node_update_tasks.wait_for_cache_task": ANALYSIS_WORKERS,
    "analysis.tasks.node_update_tasks.delete_analysis_old_node_versions": ANALYSIS_WORKERS,
    "analysis.tasks.node_update_tasks.wait_for_node": ANALYSIS_WORKERS,
    # Annotation
    "annotation.tasks.annotate_variants.delete_annotation_run": ANNOTATION_WORKERS,
    "annotation.tasks.annotate_variants.delete_annotation_run_uploaded_data": ANNOTATION_WORKERS,
    "annotation.tasks.annotate_variants.assign_range_lock_to_annotation_run": ANNOTATION_WORKERS,
    # #1654: retry clears a run's annotation rows (a big delete) then resets it to CREATED. On db_workers
    # like import_annotation_run so the bulk delete never consumes a throttled VEP slot.
    "annotation.tasks.annotate_variants.reset_annotation_run_for_retry": DB_WORKERS,
    # annotate_variants is the VEP lane (dump + VEP). Its DB upload phase is a separate task,
    # import_annotation_run, pinned to db_workers so quick bulk inserts never consume a throttled VEP
    # slot - the dispatcher runs it once a run reaches ANNOTATION_COMPLETED. See #1649.
    "annotation.tasks.annotate_variants.annotate_variants": ANNOTATION_WORKERS,
    "annotation.tasks.annotate_variants.import_annotation_run": DB_WORKERS,
    'annotation.tasks.calculate_sample_stats.calculate_vcf_stats': ANNOTATION_WORKERS,
    "annotation.tasks.cohort_sample_gene_damage_counts.CalculateCohortSampleGeneDamageCountsTask": ANNOTATION_WORKERS,
    "annotation.tasks.cohort_sample_gene_damage_counts.CohortSampleGeneDamageCountTask": ANNOTATION_WORKERS,
    "annotation.tasks.cohort_sample_gene_damage_counts.CohortSampleClassificationGeneDamageCountTask": ANNOTATION_WORKERS,
    # ClinVar import is a VCF import-processing pipeline that reuses AbstractVCFImportTaskFactory, so it
    # is routed to match the normal (genotype) VCF import: the parallel per-split data-processing lane
    # goes to WEB_WORKERS (mirrors ProcessGenotypeVCFDataTask), while the create-version and success
    # steps fall through to the default db_workers (mirrors ImportCreateVCFModelForGenotypeVCFTask /
    # ImportGenotypeVCFSuccessTask, which are likewise unrouted). It was previously all on
    # ANNOTATION_WORKERS - the VEP lane - where it serialised behind long annotate_variants runs on the
    # single-concurrency annotation worker, stalling imports for hours.
    "annotation.tasks.import_clinvar_vcf_task.ProcessClinVarVCFDataTask": WEB_WORKERS,

    # Anything that runs on data uploaded from the web should be WEB_WORKERS
    # 1. As it may be a different machine than DB workers etc.
    # 2. So that these jobs (which don't require DB access) fill up their own queue so run in parallel with db_workers
    'classification.tasks.classification_import_map_and_insert_task.ClassificationImportMapInsertTask': WEB_WORKERS,
    'snpdb.tasks.graph_generation_task.generate_graph': WEB_WORKERS,
    'upload.tasks.vcf.genotype_vcf_tasks.ProcessGenotypeVCFDataTask': WEB_WORKERS,
    'upload.tasks.import_bedfile_task.ImportBedFileTask': WEB_WORKERS,
    'upload.tasks.import_gene_coverage_task.ImportGeneCoverageTask': WEB_WORKERS,
    'upload.tasks.import_gene_list_task.ImportGeneListTask': WEB_WORKERS,
    'upload.tasks.import_patient_records_task.ImportPatientRecords': WEB_WORKERS,
    'upload.tasks.import_ped_task.ImportPedTask': WEB_WORKERS,

    # VariantID workers
    'upload.tasks.vcf.unknown_variants_task.InsertUnknownVariantsTask': VARIANT_ID_SINGLE_WORKER,
    'upload.tasks.vcf.genotype_vcf_tasks.UpdateVariantZygosityCountsTask': VARIANT_ID_SINGLE_WORKER,
    'upload.tasks.vcf.genotype_vcf_tasks.reload_vcf_task': VARIANT_ID_SINGLE_WORKER,

    # Scheduling single worker
    'analysis.tasks.analysis_update_tasks.create_and_launch_analysis_tasks': SCHEDULING_SINGLE_WORKER,
    # Periodic backlog drain (issue #346). Leases inline rather than fanning out a dispatch task
    # per analysis, so it belongs on the single worker where all leasing serialises. It's bounded
    # per tick and exits immediately when the analysis workers are already busy, so it can't
    # crowd out the interactive dispatches sharing this queue.
    'analysis.tasks.analysis_update_tasks.dispatch_analysis_backlog': SCHEDULING_SINGLE_WORKER,
    # #2667: single-authority annotation dispatcher - all run leasing/merge serialises here
    'annotation.tasks.annotation_scheduler_task.dispatch_annotation_runs': SCHEDULING_SINGLE_WORKER,
    'upload.tasks.vcf.import_vcf_step_task.schedule_pipeline_stage_steps': SCHEDULING_SINGLE_WORKER,
    'snpdb.tasks.soft_delete_tasks.remove_soft_deleted_vcfs_task': SCHEDULING_SINGLE_WORKER,

    # Partition archive
    'patients.tasks.extraction_matching_tasks.reconcile_pending_extractions': DB_WORKERS,
    'snpdb.tasks.partition_archive_tasks.perform_partition_archive': DB_WORKERS,
    "snpdb.tasks.sub_cohort_tasks.build_sub_cohort_any_sample_called_vc_task": DB_WORKERS,
}

CELERY_IMPORTS = (
    'analysis.tasks.analysis_update_tasks',
    'analysis.tasks.auto_analysis_tasks',
    'analysis.tasks.karyomapping_tasks',
    'analysis.tasks.node_update_tasks',
    'analysis.tasks.reanalysis_tasks',
    'annotation.tasks.annotate_variants',
    'annotation.tasks.annotation_scheduler_task',
    'annotation.tasks.calculate_sample_stats',
    'annotation.tasks.cohort_sample_gene_damage_counts',
    'annotation.tasks.import_clinvar_vcf_task',
    'classification.tasks.classification_import_process_variants_task',
    'classification.tasks.classification_import_task',
    'classification.tasks.classification_candidate_search_tasks',
    'classification.tasks.condition_text_automatch_task',
    'genes.tasks.gene_coverage_tasks',
    'patients.tasks.extraction_matching_tasks',
    'pedigree.models',
    'seqauto.tasks.gold_summary_tasks',
    'snpdb.models',
    'snpdb.tasks.clingen_tasks',
    'snpdb.tasks.cohort_genotype_tasks',
    'snpdb.tasks.sub_cohort_tasks',
    'snpdb.tasks.graph_generation_task',
    'snpdb.tasks.partition_archive_tasks',
    'snpdb.tasks.soft_delete_tasks',
    'snpdb.tasks.vcf_bed_file_task',
    'snpdb.tasks.vcf_zygosity_count_tasks',
    'sync.tasks.sync_tasks',
    'upload.tasks.import_bedfile_task',
    'upload.tasks.import_dragen_tso500_all_fusions_task',
    'upload.tasks.import_gene_coverage_task',
    'upload.tasks.import_gene_list_task',
    'upload.tasks.import_patient_records_task',
    'upload.tasks.import_ped_task',
    'upload.tasks.import_variant_tags_task',
    'upload.tasks.import_wiki_task',
    'upload.tasks.vcf.genotype_vcf_tasks',
    'upload.tasks.vcf.import_sql_copy_task',
    'upload.tasks.vcf.import_vcf_step_task',
    'upload.tasks.vcf.import_vcf_tasks',
    'upload.tasks.vcf.unknown_variants_task',
    'variantopedia.tasks.server_status_tasks',
)

CELERY_TASK_ALWAYS_EAGER = False  # True to execute in http server process (or Eclipse)
CELERY_RESULT_BACKEND = "redis://127.0.0.1:6379/1"

# Per-queue cap on a worker process's virtual address space (GB, RLIMIT_AS) so a runaway allocation
# raises a catchable MemoryError in Python (update_node_task perma-fails the node + reports to
# Rollbar) instead of the OS OOM-killer locking up the whole box. {} = no limit.
# A worker is capped only if EVERY queue it serves (its -Q list) is listed here, so a worker that
# also runs uncapped work (e.g. VEP bulk inserts on annotation_workers) is never throttled.
# RLIMIT_AS caps VSZ (virtual, overcounts RSS) so set generously above the worker's baseline - too
# low and the worker MemoryErrors on warmup. Tune per deployment.
CELERY_WORKER_ADDRESS_SPACE_LIMIT_GB = {
    "analysis_workers": 8,
}

# Crash safety brake: the worker records each host boot, and when two land close together (a box
# that came up then went down again - a likely crash loop) it auto-pauses the analysis + annotation
# job dispatchers, so the jobs that may have crashed the box don't immediately re-launch and crash
# it again. A single reboot the box survives starts up as normal. An admin resumes a paused
# deployment with 'manage.py jobs_control resume'. Turn this OFF on ephemeral / autoscaled hosts,
# where a fresh boot is routine rather than a crash signal.
JOBS_AUTOPAUSE_ON_REBOOT = True
JOBS_AUTOPAUSE_ON_REBOOT_UPTIME_SECS = 600  # uptime under this on worker start => treat as a reboot
JOBS_AUTOPAUSE_ON_REBOOT_WINDOW_SECS = 3600  # 2nd boot this soon after the previous => crash loop
