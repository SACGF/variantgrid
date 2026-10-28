# annotation — agent notes
Owns: AnnotationVersion + sub-versions (VariantAnnotationVersion, GeneAnnotationVersion, ClinVarVersion, HPA); the VEP
pipeline (AnnotationRangeLock → AnnotationRun → partitioned VariantAnnotation / VariantTranscriptAnnotation /
VariantGeneOverlap); AnnotSV and gene-level pipelines via AnnotationPipelineVersion; ClinVar import + ClinVarRecordCollection
XML cache; GeneAnnotation / DBNSFPGeneAnnotation; citations; cohort annotation stats. Phenotype text matching is
patients' (patients/AGENTS.md).
Start with:
- models/models_version.py (SubVersionPartition, VariantAnnotationVersion), models/models_annotation_version.py
  (AnnotationVersion), models/models_annotation_run.py (AnnotationRangeLock, AnnotationRun), models/models.py
  (VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap); ClinVar and gene/HPA sub-versions in
  models/models_clinvar.py and models/models_gene_annotation.py
- open_targets.py + vep_plugin_display.py — parsing '&'-joined VEP plugin columns for display; VariantAnnotation
  keeps same-named properties for templates
- enums/annotation_enums.py — AnnotationStatus, VariantAnnotationPipelineType, VEPPlugin / VEPCustom
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
  analyses reference (annotation/models/models_annotation_version.py:AnnotationVersion).
- Query annotation through the partition hook, never VariantAnnotation.objects directly:
  annotation/annotation_version_querysets.py:get_queryset_for_annotation_version rewrites the SQL to the version's
  partition via annotation/models/models_annotation_version.py:AnnotationVersion.sql_partition_transformer; plain joins hit the empty
  base table.
- Get "the current version" with annotation/models/models_annotation_version.py:AnnotationVersion.latest (validates, raises when none)
  or latest_or_none for degrade-gracefully paths; both filter on VAV status ACTIVE. VAVs go NEW → ACTIVE → HISTORICAL
  via annotation/models/models_version.py:VariantAnnotationVersion.promote_to_active (demotes the prior ACTIVE one and
  ERRORs its unleased incomplete runs: the dispatcher only leases runs on NEW / ACTIVE versions). Runs nothing will
  dispatch - no range lock, or a HISTORICAL version - are cleared by `clean_up_abandoned_annotation_runs`
  (annotation/abandoned_annotation_runs.py:clean_up_abandoned_annotation_runs).
- Variant pks commit in order (every variant is inserted by the one upload path, `upload/AGENTS.md`), so a new range lock
  starting above the highest locked pk can't skip a slow-to-commit insert; gaps only come from a deleted lock.
- Variants left in pk gaps between AnnotationRangeLocks (a deleted lock) are never annotated: `check_annotation_range_lock_gaps` reports them. It probes each gap row by row (`annotation/annotation_versions.py:get_range_lock_gaps_with_variants`), so it is not a deployment check: it skips the gap below the first lock (the first lock starts at the build's lowest unannotated variant) and gives each probe `--probe-timeout` seconds, listing the gaps that exceed it as unchecked.
- Saving a sub-version creates its partition tables (library/django_utils/django_partition.py:RelatedModelsPartitionModel.save)
  and bumps AnnotationVersion (annotation/models/models_version.py:SubVersionPartition.save → AnnotationVersion.new_sub_version,
  which overwrites the last one if unreferenced); batch creation goes inside SubVersionPartition.defer_new_sub_version.
- Add or gate a VEP column with a VEPColumnDef in annotation/vep_columns.py (min/max_columns_version,
  min/max_vep_version, genome_builds, pipeline_types) plus a VariantGridColumn — the mapping has no DB model.
  annotation/vep_columns.py:filter_for is the one selector the command builder, inserter and descriptions page share.
- Branch on pipeline via annotation/pipelines/__init__.py:PIPELINES (STANDARD, STRUCTURAL_VARIANT, GENE_LEVEL,
  ANNOTSV) and annotation/pipelines/base.py:AnnotationPipelineRunner (dump / annotate / import_results), not on
  pipeline_type. Which variants a type covers is annotation/annotation_pipeline_routing.py:pipeline_type_variant_q
  (pipeline_type_for_variant for one): a gene-level variant is GENE_LEVEL before any size rule (it has no size); a
  symbolic del/dup/inv under the VariantAnnotationVersion's structural_variant_min_size (pinned from
  ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE, so each version keeps its own) is STANDARD and dumped to VEP as its sequence,
  anything else symbolic is the SV pipeline, an explicit variant whose ref or alt reaches
  explicit_structural_variant_min_size is the SV pipeline (dumped as its sequence; null on VAVs made before #1708, so
  they never move one). Branch on the run's pipeline_type, never on variant.is_symbolic.
- AnnotationRun.status is recomputed from its timestamp/error fields on every save
  (annotation/models/models_annotation_run.py:AnnotationRun.get_status): set dump_start / annotation_end / upload_end, never status.
  Lifecycle: annotation_scheduler makes AnnotationRangeLocks + runs on scheduling_single_worker →
  annotation/tasks/annotation_scheduler_task.py:dispatch_annotation_runs leases (leased_by / lease_expires) and launches
  annotation/tasks/annotate_variants.py:annotate_variants on annotation_workers (dump + VEP, stops at
  ANNOTATION_COMPLETED) → import_annotation_run on db_workers bulk-inserts through
  annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:BulkVEPVCFAnnotationInserter.
- New VEP install: `manage.py create_new_variant_annotation_version`
  (annotation/annotation_versions.py:get_or_create_variant_annotation_version_from_current_vep, status NEW) → link a
  GeneAnnotationRelease (annotation/models/models_version.py:VariantAnnotationVersion.link_gene_annotation_release) →
  `gene_annotation --new-releases` → promote. Non-VEP tools: `create_new_annotation_pipeline_version`.
Gotchas:
- Call NCBI E-utilities through `library/entrez.py` (`efetch`, `esearch`, `esummary`, `epost`), keeping Bio only to parse: `Bio.Entrez.efetch` etc.
  open urllib with no timeout (`Entrez.timeout` is not read), so a stalled NCBI hangs the web worker (#2440).
- VEP's RefSeq cache names mitochondrial transcripts after the gene ('ND4.1'), never a real accession: the inserter links them to cdot's 'fake-rna-ND4' (annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:BulkVEPVCFAnnotationInserter._get_transcript_id_and_transcript_version_id), ahead of a gene-named 'ND4.1' row a pre-cdot import may have left (no genome_builds data), and stores hgvs_c as the 'NC_012920.1:m.' HGVS. tRNA/rRNA rows stay unlinked - cdot has nothing for them.
- Every VEP value in the inserter is still a string until PostgreSQL parses the COPY CSV, and gnomAD mixes notations
  ('9.1e-06', '0.25'): picking one of several values (max/min/sort) must convert with float(), or string order wins.
- A VAV must match the VEP that will run: annotation/vep_annotation.py:vep_check_command_line_version_match raises
  VEPVersionMismatchError when any data file or plugin version differs, and the annotated VCF header is checked the
  same way on import — changing a settings.ANNOTATION data path without a new VAV halts annotation.
- The first VEPVersionMismatchError latches vep_version_mismatch_reason/date on the VAV (#1679,
  annotation/models/models_version.py:VariantAnnotationVersion.latch_vep_version_mismatch); get_annotation_run_blocker then makes the
  scheduler and the dispatcher's VEP lane skip it, while the import lane keeps draining. Only
  annotation/tasks/annotate_variants.py:recheck_vep_version clears it, and get_vep_version is lru_cached per worker
  process, so the re-check clears only its own process's cache - restart annotation_workers to re-run the same VAV.
- After VEP or the bulk import, write the run with annotation/models/models_annotation_run.py:AnnotationRun.save_if_owner (conditional
  UPDATE on task_id): a lease that expired mid-VEP means another attempt owns the row and plain save() clobbers it
  (annotation/tasks/annotate_variants.py:AnnotationRunLeaseHeartbeat).
- ClinVar.variant and AnnotationRangeLock.min/max_variant are on_delete=PROTECT (annotation/models/models_clinvar.py:ClinVar,
  annotation/models/models_annotation_run.py:AnnotationRangeLock.release_variant) — clear them before deleting Variants.
- Sub-versions carry library/django_utils/data_archive_mixin.py:DataArchiveMixin; an archived partition makes
  get_variant_queryset_for_annotation_version raise DataArchivedError, so check data_archived before assuming rows.
- annotation_variantannotation partitions are per VAV and huge in production: bulk work goes through COPY / temp
  tables (annotation/backfill_columns.py:import_backfill_vcf, `annotation_backfill_columns`), never per-row ORM updates.
- Adding a VEP plugin/data file to an existing VAV (backfill instead of a new VAV): the file is a VAV pin
  (annotation/vep_config.py:vep_component_version_kwargs), so configuring it halts annotation with
  VEPVersionMismatchError and makes `create_new_variant_annotation_version` mint a full re-annotation until the
  ACTIVE VAV's new field is stamped to match. Gate the VEPColumnDef on the existing VAV's columns_version and VEP
  (min_columns_version / min_vep_version, never a columns_version bump), or `annotation_backfill_columns` refuses it.
  Its `--csq` import fills VariantAnnotation from the PICK'd entry and, for a column VariantTranscriptAnnotation
  holds, each transcript row from the entry whose Feature it is (so the VEP run keeps Feature in `--fields`).
- A backfilled_* flag on the VAV (annotation/models/models_version.py:VariantAnnotationVersion.backfilled_ptc) lets
  readers tell "not backfilled yet" from a real null: default True, the migration sets False only on VAVs that need
  the backfill, the command that fills the partition sets it True, and a ManualOperation tests for a False flag.
  The full recipe (pin stamping, gates, flag, backfill steps and limits) is claude/research/annotation.md
  "A new column on an existing VariantAnnotationVersion"; UTRAnnotator / AVI / SpliceVault (#2426) is the worked example.
- The utr_annotator_* and splicevault_* columns are transcript-level: the variant page's transcript switcher and
  evidence key autopopulate read them off the transcript rows, never VariantAnnotation, so a backfill of them has to
  fill those rows (VariantAnnotationVersion.backfilled_utr_splicevault_transcripts). A variant not yet backfilled
  shows nothing, like any other missing value.
- UTRAnnotator / AVI / SpliceVault and the AlphaGenome splicing custom VCF (#2433) are off in the package defaults. A deployment downloads (and for SpliceVault
  converts) the files, sets them in its env settings, runs `enable_vep_plugin_backfill` (stamps the pins on the
  NEW/ACTIVE versions, flags each build's newest, prints its backfill, or `--bash-script PATH` writes it as one script;
  annotation/vep_plugin_backfill.py), restarts annotation_workers, then runs the backfill - setting a file without the
  command halts the VEP lane on VEPVersionMismatchError.
- Non-commercial data (annotation/vep_config.py:NON_COMMERCIAL_VEP_CONFIG_KEYS) is cleared from VEPConfig and the pins
  when ANNOTATION_NON_COMMERCIAL_DATA_ALLOWED is False, so it never reaches the command line, the pins or the page.
  Read vep_config through VEPConfig / vep_component_version_kwargs, never genome_build.settings directly.
- The SpliceVault plugin matches the data file's transcript column against the transcript's stable_id, a versioned NM
  under --refseq, while Ensembl's file holds unversioned ENST: a RefSeq build configures the copy
  `splicevault_to_refseq` rewrites to MANE RefSeq (annotation/splicevault.py), which only covers MANE transcripts.
- annotation/models/models.py:VariantAnnotation.get_variant_annotation_q_list excludes reference (ref == alt) and '.'/'*' alt
  variants; they never get rows, so any "unannotated" count must apply it.
- annotation/models/models_annotation_run.py:AnnotationRun.get_for_variant matches by range-lock bounds and pipeline_type: a variant
  inside a lock's range with no annotation row reads as "in progress", not "missing".
- VEP `--fork` re-forks children every buffer and ships their results back over a socket, so with the pipeline's
  small `ANNOTATION_VEP_BUFFER_SIZE` it costs CPU rather than saving it; spare cores go into `annotation_workers`
  `--concurrency` instead (the comment on `ANNOTATION_VEP_FORK` in `variantgrid/settings/components/annotation_settings.py`, #2477).
- The OpenTargets VEP plugin defaults to lead_variants_only=1, which returns only the credible sets a variant leads —
  annotation/vep_annotation.py passes lead_variants_only=0 so a variant that is merely a member still gets its GWAS L2G
  score. The open_targets_* columns are '&'-joined parallel arrays zipped by
  annotation/open_targets.py:open_targets_records, so they are only ever written together.
- SVs over ANNOTATION_VEP_SV_MAX_SIZE are never dumped, so VEP's skip count never includes them: their TOO_LONG rows
  are written by annotation/vcf_files/import_vcf_annotations.py:handle_vep_skipped and, for a run with an empty dump,
  annotation/pipelines/base.py:AnnotationPipelineRunner.handle_empty_dump (#2104). Change the cap in one place:
  annotation/annotation_version_querysets.py:filter_vep_sv_max_size.
- SVGeneOverlapResolver lives in genes/gene_overlaps.py (the fusion breakpoint resolver needs it and `genes` cannot import the inserter); build it with `SVGeneOverlapResolver.for_variant_annotation_version(vav)`.
- An SV can overlap several gnomAD-SV records. annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:SVOverlapProcessor
  picks one (ANNOTATION_VEP_SV_OVERLAP_SINGLE_VALUE_METHOD) and copies its values onto the ordinary gnomad_af / gnomad_ac /
  gnomad_popmax_af / per-population columns so the analysis PopulationNode filters SVs like small variants; the
  gnomad_sv_overlap_* text fields keep them all. VEP values are still strings at that point - compare them as floats,
  as gnomAD-SV mixes '4.6e-05' and '0.006085' notation and string order puts the rarest record last.
- Reading a data-file VCF back to rewrite a column: cyvcf2 parses Float INFO through float32 ('4.6e-05' comes back as
  4.600000102072954e-05), so take the value off the raw INFO text (str(record)) to store what VEP's own pass-through would.
- Cohort stats are counted in SQL (annotation/tasks/calculate_sample_stats.py:_count_cohort_stats, #2283): rows are
  grouped by (samples_zygosity, flag bits) into an ON COMMIT DROP temp table, so it only works inside
  calculate_cohort_stats' locked transaction. Python folds the grouped rows onto the counter fields.
  A new per-row input is a bit in _BOOLEAN_FLAGS_SQL; a new counting rule goes in _flag_field_indexes / add_count.
- Trio inheritance buckets are computed apart from the cohort's (annotation/tasks/calculate_trio_stats.py, #2327) on the
  same rows and lock, keyed by the trio's members (snpdb/models/models_family.py:Trio.get_stats_filter_key). The cohort
  pass knows nothing of trios: a TrioNode miss enqueues its trio, and a base cohort recompute clears the trio rows.
  The predicates (_trio_predicates) mirror analysis/models/nodes/sources/trio_node.py.
- A sub cohort's stats are written to its parent's CGC (calculate_sample_stats.py:calculate_cohort_stats, #2311): only
  the base cohort writes per-sample rows and clears the CGC, a sub cohort writes just its aggregate rows under keys
  carrying its samples (snpdb/models/models_cohort.py:Cohort.get_stats_filter_key). Freshness needs a row per key.
- settings.ANNOTATION_GENE_ANNOTATION_VERSION_ENABLED gates whether gene_annotation_version joins the partition SQL
  (annotation/models/models_annotation_version.py:AnnotationVersion.sub_annotations_inheritance_partitioning); validate() still checks it.
- transcripts_annotation_selections.py:VariantTranscriptSelections converts the coordinate to explicit once before
  looping the other consortium's transcripts — only for alts that have an explicit form. Symbolic DEL/DUP/INV go to the
  converter as coordinates (#1571) and <CNV>/<INS> have no HGVS at all (#1574), so both are passed through symbolic.
Tests:
- annotation/fake_data.py:get_fake_annotation_version builds a full valid AnnotationVersion (ACTIVE VAV at
  columns_version 2, fake gene release, ontology, ClinVar, HPA); refuses to run outside UNIT_TEST.
- variantgrid/test_runner.py:VariantGridTestRunner seeds get_fake_annotation_version for GRCh37/GRCh38 into the
  main test database before Django clones it for the workers, so the per-class calls are lookups. A test that
  needs its own ACTIVE VariantAnnotationVersion (one_active_vav_per_build) calls
  annotation/fake_data.py:retire_seeded_annotation_version(genome_build) first.
- annotation/fake_data.py:create_fake_variants loads the fixture VCF's variants (wraps
  snpdb/tests/utils/vcf_testing_utils.py:slowly_create_loci_and_variants_for_vcf); create_fake_variant_annotation and
  create_fake_clinvar_data add rows. genes/fake_data.py:create_fake_transcript_version (and the
  create_gata2_transcript_version / create_pten_transcript_version presets) make genes/transcripts with exon structure.
- An ORM insert into VariantAnnotation lands in the base table, which analysis nodes never read - they join the
  version's partition. Rows a node must see are written inside `temporary_db_table` (fake_data.py:FakeVariants).
- @override_settings(**get_fake_annotation_settings_dict(columns_version=N)) pins settings.ANNOTATION to the
  tests/test_data/test_columns_versionN_<build>.vep_annotated.vcf fixtures; regenerate them with `vep_run --test`
  (needs a real VEP install).
- tests/test_urls.py is the URLTestCase; library/django_utils/unittest_utils.py:URLTestCase sets
  ANNOTATION_CACHED_WEB_RESOURCES=[] so CachedWebResource loaders never fire. tests/test_annotation_vcf.py drives the
  real BulkVEPVCFAnnotationInserter over every columns_version — the slow one.
Deep reference: __annotation_readme.md · claude/research/annotation.md · claude/maps/models.md#annotation
