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
BROKER_URL = get_secret("CELERY.broker_url")

CELERY_DEFAULT_QUEUE = 'db_workers'
CELERY_QUEUES = (
    Queue('analysis_workers', Exchange('analysis_workers'), routing_key='analysis_workers'),
    Queue('annotation_workers', Exchange('annotation_workers'), routing_key='annotation_workers'),
    # db_workers - only use the database, ie don't write to web server filesystem (so it could be moved to diff machine)
    Queue('db_workers', Exchange('db_workers'), routing_key='db_workers'),
    Queue('web_workers', Exchange('web_workers'), routing_key='web_workers'),
    Queue('scheduling_single_worker', Exchange('scheduling_single_worker'), routing_key='scheduling_single_worker'),
    Queue('variant_id_single_worker', Exchange('variant_id_single_worker'), routing_key='variant_id_single_worker'),
    Queue('seqauto_single_worker', Exchange('seqauto_single_worker'), routing_key='seqauto_single_worker'),
)

ANALYSIS_WORKERS = {'queue': 'analysis_workers', 'routing_key': 'analysis_workers'}
ANNOTATION_WORKERS = {'queue': 'annotation_workers', 'routing_key': 'annotation_workers'}
DB_WORKERS = {'queue': 'db_workers', 'routing_key': 'db_workers'}  # Default Queue
# Require reading from filesystem (upload directory or images & reports)
WEB_WORKERS = {'queue': 'web_workers', 'routing_key': 'web_workers'}

# This is used for inserting new Variants/Loci and updating Redis with hashes
# 1 worker to avoid race conditions so we end up with 1 and only 1 variant for chrom/post/ref/alt
VARIANT_ID_SINGLE_WORKER = {'queue': 'variant_id_single_worker', 'routing_key': 'variant_id_single_worker'}

# 1 worker, use this to schedule tasks and avoid race conditions
SCHEDULING_SINGLE_WORKER = {'queue': 'scheduling_single_worker', 'routing_key': 'scheduling_single_worker'}

SEQAUTO_SINGLE_WORKERS = {'queue': 'seqauto_single_worker', 'routing_key': 'seqauto_single_worker'}

CELERY_WORKER_NAMES = ['annotation_workers', 'db_workers', 'web_workers',
                       'scheduling_single_worker', 'variant_id_single_worker']
CELERY_ANALYSIS_WORKER_NAMES = ['analysis_workers']
CELERY_SEQAUTO_WORKER_NAMES = ['seqauto_single_worker']

CELERY_ROUTES = {
    # Analysis
    'analysis.tasks.karyomapping_tasks.create_genome_karyomapping_for_trio': ANALYSIS_WORKERS,
    "analysis.tasks.node_update_tasks.dummy_task": ANALYSIS_WORKERS,
    "analysis.tasks.node_update_tasks.update_node_task": ANALYSIS_WORKERS,
    "analysis.tasks.node_update_tasks.node_cache_task": ANALYSIS_WORKERS,
    "analysis.tasks.node_update_tasks.wait_for_cache_task": ANALYSIS_WORKERS,
    "analysis.tasks.node_update_tasks.delete_old_node_versions": ANALYSIS_WORKERS,
    "analysis.tasks.node_update_tasks.wait_for_node": ANALYSIS_WORKERS,
    # Annotation
    "annotation.tasks.annotate_variants.delete_annotation_run": ANNOTATION_WORKERS,
    "annotation.tasks.annotate_variants.delete_annotation_run_uploaded_data": ANNOTATION_WORKERS,
    "annotation.tasks.annotate_variants.assign_range_lock_to_annotation_run": ANNOTATION_WORKERS,
    "annotation.tasks.annotate_variants.annotate_variants": ANNOTATION_WORKERS,
    'annotation.tasks.calculate_sample_stats.calculate_vcf_stats': ANNOTATION_WORKERS,
    "annotation.tasks.cohort_sample_gene_damage_counts.CalculateCohortSampleGeneDamageCountsTask": ANNOTATION_WORKERS,
    "annotation.tasks.cohort_sample_gene_damage_counts.CohortSampleGeneDamageCountTask": ANNOTATION_WORKERS,
    "annotation.tasks.cohort_sample_gene_damage_counts.CohortSampleClassificationGeneDamageCountTask": ANNOTATION_WORKERS,
    "annotation.tasks.import_clinvar_vcf_task.ImportCreateVersionForClinVarVCFTask": ANNOTATION_WORKERS,
    "annotation.tasks.import_clinvar_vcf_task.ProcessClinVarVCFDataTask": ANNOTATION_WORKERS,
    "annotation.tasks.import_clinvar_vcf_task.ImportClinVarSuccessTask": ANNOTATION_WORKERS,

    # Anything that runs on data uploaded from the web should be WEB_WORKERS
    # 1. As it may be a different machine than DB workers etc.
    # 2. So that these jobs (which don't require DB access) fill up their own queue so run in parallel with db_workers
    'upload.tasks.vcf.import_vcf_task.ImportVCFExtractUnknownVariantsAndSplitFileTask': WEB_WORKERS,
    'upload.tasks.vcf.genotype_vcf_tasks.ImportCreateProjectForGenotypeVCFTask': WEB_WORKERS,
    'upload.tasks.vcf.genotype_vcf_tasks.ProcessGenotypeVCFDataTask': WEB_WORKERS,
    'upload.tasks.import_bedfile_task.ImportBedFileTask': WEB_WORKERS,
    'upload.tasks.import_expression_task.ImportExpressionTask': WEB_WORKERS,
    'upload.tasks.import_gene_list_task.ImportGeneCoverageTask': WEB_WORKERS,
    'upload.tasks.import_gene_list_task.ImportGeneListTask': WEB_WORKERS,
    'upload.tasks.import_patient_records_task.ImportPatientRecords': WEB_WORKERS,
    'upload.tasks.import_ped_task.ImportPedTask': WEB_WORKERS,
    'generate_graph': WEB_WORKERS,

    # VariantID workers
    'upload.tasks.vcf.unknown_variants_task.InsertUnknownVariantsTask': VARIANT_ID_SINGLE_WORKER,
    'upload.tasks.load_variants_hash_in_redis_task.load_variants_hash_in_redis': VARIANT_ID_SINGLE_WORKER,
    'upload.tasks.vcf.genotype_vcf_tasks.UpdateVariantZygosityCountsTask': VARIANT_ID_SINGLE_WORKER,
    'upload.tasks.vcf.genotype_vcf_tasks.reload_vcf_task': VARIANT_ID_SINGLE_WORKER,
    'snpdb.tasks.validation_task.validate_variant_data': VARIANT_ID_SINGLE_WORKER,

    # Scheduling single worker
    'upload.tasks.annotation_scheduler_task.annotation_scheduler': SCHEDULING_SINGLE_WORKER,
    'upload.tasks.vcf.import_vcf_step_task.schedule_pipeline_stage_steps': SCHEDULING_SINGLE_WORKER,
    'snpdb.tasks.soft_delete_tasks.remove_soft_deleted_vcfs_task': SCHEDULING_SINGLE_WORKER,

    # SeqAuto
    'seqauto.tasks.scan_run_jobs.scan_run_jobs': SEQAUTO_SINGLE_WORKERS,
    'sapath.tasks.import_helix_task.sapath_helix_load_if_changed': SEQAUTO_SINGLE_WORKERS,
}

CELERY_IMPORTS = (
    'analysis.tasks.karyomapping_tasks',
    'analysis.tasks.node_update_tasks',
    'annotation.tasks.annotate_variants',
    'annotation.tasks.calculate_sample_stats',
    'annotation.tasks.cohort_sample_gene_damage_counts',
    'annotation.tasks.import_clinvar_vcf_task',
    'genes.tasks.gene_coverage_tasks',
    'pedigree.models',
    'seqauto.tasks.gold_summary_tasks',
    'seqauto.tasks.scan_run_jobs',
    'snpdb.models',
    'snpdb.tasks.clingen_tasks',
    'snpdb.tasks.cohort_genotype_tasks',
    'snpdb.tasks.graph_generation_task',
    'snpdb.tasks.soft_delete_tasks',
    'snpdb.tasks.vcf_bed_file_task',
    'snpdb.tasks.validation_task',
    'upload.tasks.import_bedfile_task',
    'upload.tasks.import_expression_task',
    'upload.tasks.import_gene_coverage_task',
    'upload.tasks.import_gene_list_task',
    'upload.tasks.import_patient_records_task',
    'upload.tasks.import_ped_task',
    'upload.tasks.import_variant_tags_task',
    'upload.tasks.load_variants_hash_in_redis_task',
    'upload.tasks.vcf.genotype_vcf_tasks',
    'upload.tasks.vcf.import_sql_copy_task',
    'upload.tasks.vcf.import_vcf_step_task',
    'upload.tasks.vcf.import_vcf_tasks',
    'upload.tasks.vcf.unknown_variants_task',
    'classification.tasks.classification_import_process_variants_task',
    'classification.tasks.classification_import_task',
    'sync.tasks.sync_tasks'
)

CELERY_ALWAYS_EAGER = False  # True to execute in http server process (or Eclipse)
CELERY_RESULT_BACKEND = BROKER_URL
