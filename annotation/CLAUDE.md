# annotation — agent notes
Owns: AnnotationVersion + sub-versions (VariantAnnotationVersion, GeneAnnotationVersion, ClinVarVersion, HPA); the VEP
pipeline (AnnotationRangeLock → AnnotationRun → partitioned VariantAnnotation / VariantTranscriptAnnotation /
VariantGeneOverlap); AnnotSV and gene-level pipelines via AnnotationPipelineVersion; ClinVar import + ClinVarRecordCollection
XML cache; GeneAnnotation / DBNSFPGeneAnnotation; citations; cohort annotation stats; phenotype text matching.
Start with:
- models/models.py — every version model, AnnotationRun, VariantAnnotation (~3,000 lines: outline with grep)
- models/models_enums.py — AnnotationStatus, VariantAnnotationPipelineType, VEPPlugin / VEPCustom
- annotation_version_querysets.py — how to query one version's partition
- vep_columns.py — VEPColumnDef registry: VEP output field → VariantGridColumn, gated by columns/VEP version
- vep_annotation.py + vep_config.py — the VEP command line and settings.ANNOTATION[build] as VEPConfig
- tasks/annotation_scheduler_task.py + tasks/annotate_variants.py — range locks, dispatch, lease, VEP lane, import lane
- pipelines/__init__.py — PIPELINES registry, one runner per pipeline type
Patterns here:
- Keep the three "versions" apart. settings.ANNOTATION_VEP_VERSION is the VEP binary. VariantAnnotationVersion (VAV)
  is a DB snapshot of that VEP plus every data file/plugin it ran with, one per build; its columns_version selects
  which VEPColumnDefs apply (annotation/vep_config.py:VEPConfig reads genome_build.settings["columns_version"]).
  AnnotationVersion bundles a VAV with the gene/ClinVar/HPA/ontology sub-versions and is what querysets and
  analyses reference (annotation/models/models.py:AnnotationVersion).
- Query annotation through the partition hook, never VariantAnnotation.objects directly:
  annotation/annotation_version_querysets.py:get_queryset_for_annotation_version rewrites the SQL to the version's
  partition via annotation/models/models.py:AnnotationVersion.sql_partition_transformer; plain joins hit the empty
  base table.
- Get "the current version" with annotation/models/models.py:AnnotationVersion.latest (validates, raises when none)
  or latest_or_none for degrade-gracefully paths; both filter on VAV status ACTIVE. VAVs go NEW → ACTIVE → HISTORICAL
  via annotation/models/models.py:VariantAnnotationVersion.promote_to_active (demotes the prior ACTIVE one).
- Saving a sub-version creates its partition tables (library/django_utils/django_partition.py:RelatedModelsPartitionModel.save)
  and bumps AnnotationVersion (annotation/models/models.py:SubVersionPartition.save → AnnotationVersion.new_sub_version,
  which overwrites the last one if unreferenced); batch creation goes inside SubVersionPartition.defer_new_sub_version.
- Add or gate a VEP column with a VEPColumnDef in annotation/vep_columns.py (min/max_columns_version,
  min/max_vep_version, genome_builds, pipeline_types) plus a VariantGridColumn — the mapping has no DB model.
  annotation/vep_columns.py:filter_for is the one selector the command builder, inserter and descriptions page share.
- Branch on pipeline via annotation/pipelines/__init__.py:PIPELINES (STANDARD, STRUCTURAL_VARIANT, GENE_LEVEL,
  ANNOTSV) and annotation/pipelines/base.py:AnnotationPipelineRunner (dump / annotate / import_results), not on
  pipeline_type. Which variants a type covers is annotation/annotation_version_querysets.py:pipeline_type_variant_q:
  symbolic (svlen) variants are the SV pipeline, gene-level variants belong to neither VEP pipeline.
- AnnotationRun.status is recomputed from its timestamp/error fields on every save
  (annotation/models/models.py:AnnotationRun.get_status): set dump_start / annotation_end / upload_end, never status.
  Lifecycle: annotation_scheduler makes AnnotationRangeLocks + runs on scheduling_single_worker →
  annotation/tasks/annotation_scheduler_task.py:dispatch_annotation_runs leases (leased_by / lease_expires) and launches
  annotation/tasks/annotate_variants.py:annotate_variants on annotation_workers (dump + VEP, stops at
  ANNOTATION_COMPLETED) → import_annotation_run on db_workers bulk-inserts through
  annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:BulkVEPVCFAnnotationInserter.
- New VEP install: `manage.py create_new_variant_annotation_version`
  (annotation/annotation_versions.py:get_or_create_variant_annotation_version_from_current_vep, status NEW) → link a
  GeneAnnotationRelease (annotation/models/models.py:VariantAnnotationVersion.link_gene_annotation_release) →
  `gene_annotation --new-releases` → promote. Non-VEP tools: `create_new_annotation_pipeline_version`.
Gotchas:
- A VAV must match the VEP that will run: annotation/vep_annotation.py:vep_check_command_line_version_match raises
  VEPVersionMismatchError when any data file or plugin version differs, and the annotated VCF header is checked the
  same way on import — changing a settings.ANNOTATION data path without a new VAV halts annotation.
- After VEP or the bulk import, write the run with annotation/models/models.py:AnnotationRun.save_if_owner (conditional
  UPDATE on task_id): a lease that expired mid-VEP means another attempt owns the row and plain save() clobbers it
  (annotation/tasks/annotate_variants.py:AnnotationRunLeaseHeartbeat).
- ClinVar.variant and AnnotationRangeLock.min/max_variant are on_delete=PROTECT (annotation/models/models.py:ClinVar,
  annotation/models/models.py:AnnotationRangeLock.release_variant) — clear them before deleting Variants.
- Sub-versions carry library/django_utils/data_archive_mixin.py:DataArchiveMixin; an archived partition makes
  get_variant_queryset_for_annotation_version raise DataArchivedError, so check data_archived before assuming rows.
- annotation_variantannotation partitions are per VAV and huge in production: bulk work goes through COPY / temp
  tables (annotation/backfill_columns.py:import_backfill_vcf, `annotation_backfill_columns`), never per-row ORM updates.
- annotation/models/models.py:VariantAnnotation.VARIANT_ANNOTATION_Q excludes reference (ref == alt) and '.'/'*' alt
  variants; they never get rows, so any "unannotated" count must apply it.
- annotation/models/models.py:AnnotationRun.get_for_variant matches by range-lock bounds and pipeline_type: a variant
  inside a lock's range with no annotation row reads as "in progress", not "missing".
- settings.ANNOTATION_GENE_ANNOTATION_VERSION_ENABLED gates whether gene_annotation_version joins the partition SQL
  (annotation/models/models.py:AnnotationVersion.sub_annotations_inheritance_partitioning); validate() still checks it.
Tests:
- annotation/fake_annotation.py:get_fake_annotation_version builds a full valid AnnotationVersion (ACTIVE VAV at
  columns_version 2, fake gene release, ontology, ClinVar, HPA); refuses to run outside UNIT_TEST.
- annotation/fake_annotation.py:create_fake_variants loads the fixture VCF's variants (wraps
  snpdb/tests/utils/vcf_testing_utils.py:slowly_create_loci_and_variants_for_vcf); create_fake_variant_annotation and
  create_fake_clinvar_data add rows. annotation/tests/test_data_fake_genes.py:create_fake_transcript_version (and the
  create_gata2_transcript_version / create_pten_transcript_version presets) make genes/transcripts with exon structure.
- @override_settings(**get_fake_annotation_settings_dict(columns_version=N)) pins settings.ANNOTATION to the
  tests/test_data/test_columns_versionN_<build>.vep_annotated.vcf fixtures; regenerate them with `vep_run --test`
  (needs a real VEP install).
- tests/test_urls.py is the URLTestCase; library/django_utils/unittest_utils.py:URLTestCase sets
  ANNOTATION_CACHED_WEB_RESOURCES=[] so CachedWebResource loaders never fire. tests/test_annotation_vcf.py drives the
  real BulkVEPVCFAnnotationInserter over every columns_version — the slow one.
Deep reference: __annotation_readme.md · claude/research/annotation.md · claude/maps/models.md#annotation
