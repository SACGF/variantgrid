# annotation — research notes

Verified against 7c4408c62 on 2026-09-06

The annotation app takes the Variants snpdb has stored and attaches what the world knows about them: VEP's per-transcript
consequences and the plugin/custom scores that ride along with it, ClinVar's summary per variant, per-gene annotation
(ontology terms, GenCC gene-disease strength, dbNSFP gene scores, gnomAD constraint), Human Protein Atlas tissue data,
and citations. Every one of these is versioned, and `annotation/models/models.py:AnnotationVersion` is the bundle per
build that an analysis or a variant page pins itself to. This document is the story behind `annotation/CLAUDE.md`:
how a version comes to exist, how a variant gets its rows, why the tables are partitioned, and what has gone wrong
before. Model fields, URLs, commands, tasks and signals are in the generated maps
([models](../maps/models.md#annotation), [urls](../maps/urls.md#annotation), [commands](../maps/commands.md),
[tasks](../maps/tasks.md#annotation), [settings](../maps/settings.md)); `claude/domain.md` has the vocabulary.

## Flows

### A new VariantAnnotationVersion

A `annotation/models/models.py:VariantAnnotationVersion` (VAV) is a database snapshot of one VEP install for one build:
the VEP code version, the cache version, and a pin for every data file and plugin the command line will reference
(gnomAD, dbNSFP, SpliceAI, COSMIC, MaveDB, Open Targets ...), plus `columns_version`, which selects the set of columns
the pipeline writes. `annotation/vep_annotation.py:get_vep_version` runs `vep --help`-style introspection and
`annotation/vep_annotation.py:vep_dict_to_variant_annotation_version_kwargs` turns that plus
`annotation/vep_config.py:VEPConfig` (the `settings.ANNOTATION[build]` block) into the kwargs;
`annotation/annotation_versions.py:get_or_create_variant_annotation_version_from_current_vep` then `get_or_create`s the
row with `status=NEW`, which is what `manage.py create_new_variant_annotation_version` does per build. Nothing is
re-annotated yet: NEW is a version being built.

The VAV must know which gene set VEP used, because the import lane links every annotation row to a
`TranscriptVersion` by FK. `annotation/gene_annotation_release_matching.py:VepGeneSetVersions` reads the release
token out of the VAV's `refseq` / `genebuild` / `vep_cache` pins (eg `RS_2025_08`), and
`annotation/models/models.py:VariantAnnotationVersion.link_gene_annotation_release` looks for an existing
`GeneAnnotationRelease` whose cdot import URL names that release; if there is none, the genes app's
`import_cdot_gene_annotation_release` command downloads the matching cdot per-GFF file and creates it (see
`claude/research/genes.md`). A VAV with no release is a VAV whose every run will fail at pipeline time, which is what
`annotation/models/models.py:VariantAnnotationVersion.get_annotation_run_blocker` reports on the runs page. Gene
annotation is then built for that release (`manage.py gene_annotation --new-releases`, below), and only then can
`annotation/models/models.py:VariantAnnotationVersion.promote_to_active` run: it refuses while
`get_gene_annotation_promote_blocker` says the AnnotationVersion would be inconsistent, then demotes the prior ACTIVE
VAV for the build to HISTORICAL inside a `select_for_update`. A partial unique constraint (`one_active_vav_per_build`)
makes "the VAV" a database fact rather than a convention.

Saving any sub-version creates its partition tables (`annotation/models/models.py:SubVersionPartition.save` via
`library/django_utils/django_partition.py:RelatedModelsPartitionModel`) and calls
`annotation/models/models.py:AnnotationVersion.new_sub_version`, which assembles a new AnnotationVersion from the latest
sub-version of each kind for the build - overwriting the previous AnnotationVersion in place if nothing references it
yet, because a deploy that imports ClinVar, HPA and a new VAV in one go should not leave three half-built versions
behind. `SubVersionPartition.defer_new_sub_version` holds the bump back while a batch of sub-versions is written.
`AnnotationVersion.validate` is the consistency check the rest of the system relies on: every sub-version present, the
GeneAnnotationVersion built for the same GeneAnnotationRelease as the VAV, and the same OntologyVersion.

### Range locks, runs and the two lanes

Annotation is scheduled by variant pk range. `annotation/tasks/annotation_scheduler_task.py:annotation_scheduler`
(beat, and kicked by `upload/tasks/vcf/import_vcf_tasks.py:CheckStartAnnotationTask` after a VCF's variants are
inserted) runs on `scheduling_single_worker` under a cache lock and, for each annotated build's latest VAV of the
requested status, loops `annotation/tasks/annotation_scheduler_task.py:_handle_variant_annotation_version` until
`annotation/annotation_versions.py:get_annotation_range_lock_and_unannotated_count` finds no more unannotated variants.
That function walks pk blocks above the highest existing lock (never a sort-and-limit over the whole table, which was
the original design and stopped scaling) and returns an `annotation/models/models.py:AnnotationRangeLock` of at least
`ANNOTATION_VEP_BATCH_MIN` and at most `ANNOTATION_VEP_BATCH_MAX` variants. "Unannotated" is defined once, in
`annotation/annotation_version_querysets.py:get_variants_qs_for_annotation`: variants in the build's contigs passing
`annotation/models/models.py:VariantAnnotation.VARIANT_ANNOTATION_Q` (no reference-alt, no `.`/`*` alt) with no
VariantAnnotation row in this version's partition, restricted to a pipeline type by
`annotation/annotation_version_querysets.py:pipeline_type_variant_q`.

Each lock gets one `annotation/models/models.py:AnnotationRun` per enabled pipeline type
(`annotation/tasks/annotation_scheduler_task.py:_handle_range_lock`). The types are the registry in
`annotation/pipelines/__init__.py:PIPELINES`: STANDARD (VEP over short variants), STRUCTURAL_VARIANT (VEP over symbolic
variants), GENE_LEVEL (fusions, computed locally by `annotation/gene_level_annotation.py:annotate_gene_level_run`) and
ANNOTSV, which depends on the SV run having FINISHED and updates the rows it wrote. The scheduler only creates pending
state. `annotation/tasks/annotation_scheduler_task.py:dispatch_annotation_runs` is the single authority that launches:
it sweeps the latest NEW and ACTIVE VAV per build, reclaims stalled runs, fills the count lane
(`annotation/tasks/annotation_scheduler_task.py:count_annotation_runs` pre-counts each run so an empty one finishes
without a dump), merges consecutive pending locks into bigger batches once workers are busy
(`annotation/annotation_versions.py:merge_pending_range_locks`, the inverse of
`annotation/tasks/annotation_scheduler_task.py:subdivide_annotation_range_lock`), and leases up to the free capacity
with `annotation/tasks/annotation_scheduler_task.py:_lease_and_launch_run`. Leases (`leased_by`, `lease_expires`) are
written only on `scheduling_single_worker`, so they serialise without row-lock gymnastics.

A run has two stages on two queues. `annotation/tasks/annotate_variants.py:annotate_variants` (`annotation_workers`)
grabs the `task_id` lock with a conditional UPDATE, checks the run blocker, and inside an
`annotation/tasks/annotate_variants.py:AnnotationRunLeaseHeartbeat` calls the runner's `annotate`:
`annotation/pipelines/vep.py:VEPRunner.annotate` dumps the range with `annotation/annotation_run_files.py:write_qs_to_vcf`
(the VCF carries `variant_id` in INFO so the import never has to re-resolve coordinates), checks the command-line VEP
against the VAV's pins (`annotation/vep_annotation.py:vep_check_command_line_version_match`), runs
`annotation/vep_annotation.py:get_vep_command` under a memory cap, and stops at ANNOTATION_COMPLETED. The dispatcher then
re-picks the run in its resume lane and launches `annotation/tasks/annotate_variants.py:import_annotation_run` on
`db_workers`, which checks the annotated file's header the same way, streams it through
`annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:BulkVEPVCFAnnotationInserter`, clears any stale error and sends
`annotation_run_complete_signal`. Status is never set directly: `annotation/models/models.py:AnnotationRun.get_status`
derives it from the timestamp and error fields on every save, and `annotation/models/models.py:AnnotationRun.save_if_owner`
is the write used after any stage during which the lease could have expired.

### Import: VEP VCF to partition rows

`BulkVEPVCFAnnotationInserter.__init__` reads the CSQ header, selects the applicable
`annotation/vep_columns.py:VEPColumnDef` entries with `annotation/vep_columns.py:filter_for` (build, pipeline type,
columns_version, VEP version, COSMIC release, gnomAD 4 minor version, and whether the data file is actually configured),
and validates that every expected source field is in the header. `process_entry` splits the CSQ into one record per
transcript, links each to a Transcript / TranscriptVersion pk by accession (exact version only - the inserter builds its
`HGVSMatcher` with `allow_alternative_transcript_version=False`), picks the record VEP flagged with `--flag_pick` as the
representative `VariantAnnotation`, writes every record as a `VariantTranscriptAnnotation`, and one `VariantGeneOverlap`
per overlapping gene. Calculated columns (SpliceAI max, prediction counts, MaxEntScan diff, PTC/NMD distances #579,
COSMIC id merge, local c.HGVS and g.HGVS) are added in `add_calculated_variant_annotation_columns` and
`add_calculated_transcript_columns`. `bulk_insert` writes CSVs and `COPY`s them straight into the VAV's partition
tables; under `UNIT_TEST` it falls back to `bulk_create` into the base table, which is why tests never exercise the
partition path.

Reading goes the other way through one hook. Partitions are old-style inheritance children, so a plain ORM join to
`annotation_variantannotation` hits the empty parent; `annotation/annotation_version_querysets.py:get_queryset_for_annotation_version`
attaches `annotation/models/models.py:AnnotationVersion.sql_partition_transformer` to the queryset via
`library/django_utils/django_queryset_sql_transformer.py:get_queryset_with_transformer_hook`, and the compiled SQL has
its table names rewritten to the version's partitions. Django's `FilteredRelation` was considered and rejected because
it cannot follow nested relations (`variantannotation__gene__geneannotation`).

### ClinVar

ClinVar arrives two ways. The bulk summary is the ClinVar VCF, imported through the upload app's step pipeline
(`annotation/import_task_factories.py` registers the steps): `annotation/tasks/import_clinvar_vcf_task.py:ImportCreateVersionForClinVarVCFTask`
makes a `annotation/models/models.py:ClinVarVersion` keyed on the file's sha256 (a re-upload of the same file wipes and
refills the same version's partition), unknown variants are inserted by the normal VCF machinery, then
`annotation/vcf_files/import_clinvar_vcf.py:BulkClinVarInserter` maps INFO fields onto `annotation/models/models.py:ClinVar`
rows, including the oncogenic (ONC*) and somatic clinical impact (SCI*) fields added in #1791. Only staff may upload one
(`annotation/vcf_files/import_clinvar_vcf.py:check_can_import_clinvar`). The per-submission detail is fetched lazily
from ClinVar's XML: `annotation/clinvar_fetch_request.py:ClinVarFetchRequest.fetch` caches a
`annotation/models/models.py:ClinVarRecordCollection` per variation id, refreshing when older than
`CLINVAR_RECORD_CACHE_DAYS` or when `annotation/clinvar_xml_parser_via_vcv.py:ClinVarXmlParserViaVCV` bumps its
`PARSER_VERSION`, holding a `select_for_update` so two page loads do not both hit the API.

### Gene annotation

`annotation/models/models.py:GeneAnnotation` is a per-gene, per-`annotation/models/models.py:GeneAnnotationVersion` row
of pipe-joined OMIM/HPO/MONDO terms, GenCC strength flags, a dbNSFP gene row and gnomAD LoF o/e. It is built by
`annotation/management/commands/gene_annotation.py:Command`, which must be re-run whenever the GeneAnnotationRelease or
the OntologyVersion changes because the version is unique on both: `_populate_gene_annotation_version` matches every
symbol into the release with `genes/gene_matching.py:ReleaseGeneMatcher`, refreshes PanelApp Australia in one crawl,
then walks every HGNC term through an `ontology/ontology_traversal.py:get_ontology_traverser` (in-memory graph optional)
and groups the results by gene id before a bulk write. A version that matches nothing raises rather than leaving an
empty version live, since an empty GeneAnnotation partition makes every OMIM/HPO filter silently match nothing.
`--missing` and `--new-releases` are the deploy-time modes; `--force` rebuilds.

## Why it is shaped this way

Three "versions" are kept apart on purpose (`annotation/CLAUDE.md` says how to tell them apart). The VAV is what VEP
ran with, and its partition is the unit of retention: old analyses keep reading their old partition after an upgrade,
and dropping a HISTORICAL VAV drops one table rather than deleting millions of rows. `columns_version` exists so that a
VEP upgrade that changes nothing in the schema is a new VAV without a migration, and a schema change (new plugin
columns) is a new `VEPColumnDef` gated by `min_columns_version`; `annotation/vep_columns.py:VEP_COLUMNS` is code, not a
table, so the mapping is versioned with the code that reads it. `annotation/models/models.py:AnnotationPipelineVersion`
(#720) is deliberately not the VAV: pinning AnnotSV on the VEP version would make rolling its data bundle cost a full
re-annotation, so non-VEP tools get their own NEW/ACTIVE/HISTORICAL lifecycle and their runs hang off the locks the VEP
run already created - "this lock has no run at the ACTIVE version" is the same query as "this lock has no run of this
type at all", so backfill and upgrade are one code path.

The scheduler/dispatcher split (#2667) mirrors the analysis app's node tasks: creating state and launching it are
different jobs, and only the launcher needs to know about capacity. Latency comes first - while workers are free,
pending work launches as-is - and merging only kicks in when a backlog forms, so a single uploaded VCF is annotated
promptly while a bulk import ends up in efficient batches. The VEP lane and the import lane are separate tasks on
separate queues (#1649) because VEP throttling is about memory (`ANNOTATION_VEP_MEMORY_LIMIT_GB`, #1710) and a quick
bulk insert should never occupy a VEP slot. Leases with a heartbeat (#1658) exist because VEP can run for hours on
structural variants; a dead worker's run is reclaimed after `ANNOTATION_RUN_LEASE_SECONDS`, and every write after VEP
is a conditional UPDATE on `task_id` so the losing attempt cannot clobber the run its successor now owns.
`attempt_count` is bumped at execution rather than dispatch so a run that merely sat in a starved queue never burns its
`ANNOTATION_MAX_RUN_ATTEMPTS`.

External annotation (#1568, `annotation/external_annotation.py`) exists because production VEP runs were too slow on
the web host: runs on a NEW VAV can be dumped as self-describing gzipped VCFs with a sidecar of the VAV identity and
the exact VEP command templated for another machine, and imported back into the normal import lane. The range endpoints
are verified against local Variant pks on the way back (`annotation/external_annotation.py:verify_annotated_vcf_variant_ids`)
because the dump's `variant_id` INFO is the only link between the external file and the database.

## History

The app began as a flat "annotate everything with VEP into one table" design; partitioning per VAV came with the need to
keep old analyses reproducible. `ColumnVEPField` was a database table mapping VEP fields to columns until migration
`annotation/migrations/0129_delete_columnvepfield.py` replaced it with the `VEPColumnDef` registry in
`annotation/vep_columns.py`, after data-driven column gating proved impossible to test or diff. `columns_version` 2
brought dbNSFP rankscores, 3 gnomAD v4 and the `_xy` fields, 4 raw dbNSFP 5.3 scores with per-transcript resolution
(`annotation/refseq_ensembl_resolver.py`), denovo-db and the VEP 116 plugins (#1652). The VAV `active` boolean became
the NEW/ACTIVE/HISTORICAL `status` in `annotation/migrations/0143_variantannotationversion_status.py` (#577) so a version
could be built without being live. AnnotSV was split out of the SV VEP run into its own pipeline type and then given
`AnnotationPipelineVersion` (#720, migration `annotation/migrations/0168_annotation_pipeline_version_lifecycle.py`).
Truncated VEP output was caught by counting VEP's `--skipped_variants_file` against the dump (#1701); SV conservation
scoring moved from a VEP `--custom` bigWig to pyBigWig (#1657); gene-level (fusion) variants got their own local
pipeline (#1506). PyHGVS was removed in favour of biocommons throughout (#1678). Column backfill from an annotated VCF
without re-running the pipeline (`annotation/backfill_columns.py:import_backfill_vcf`, #1675) and Open Targets L2G
scores (#1822) are the most recent additions.

## Traps

The pins on a VAV must match the VEP that will run and the VCF that comes back: `vep_check_command_line_version_match`
and `annotation/vep_annotation.py:vep_check_annotated_file_version_match` both raise `VEPVersionMismatchError`, so
changing a `settings.ANNOTATION` data path without creating a new VAV halts every run. `gencode_subset` and `distance`
are settings snapshots, not VEP header facts, and are excluded from the comparison
(`annotation/vep_annotation.py:_vep_check_version_match`).

`annotation/models/models.py:AnnotationRun.get_for_variant` answers "is this variant in progress" by range-lock bounds
and pipeline type, so a variant inside a lock with no row reads as pending, not missing; an unannotated count must apply
`VARIANT_ANNOTATION_Q` or it will count reference variants that never get rows. Range locks and ClinVar rows PROTECT
their variants (`annotation/models/models.py:AnnotationRangeLock.release_variant` is the way to move a lock's endpoint
off a variant you need to delete).

Merging is only safe while every run on a lock is still CREATED and un-leased
(`annotation/annotation_versions.py:_range_lock_is_dispatchable`); an empty-finished run is reopened when its lock grows
because it owns no rows (`annotation/models/models.py:AnnotationRun.reopen_to_created`). A retry must reuse the run in
place (`annotation/models/models.py:AnnotationRun.reset_for_retry`, #1654) - deleting and recreating left rangeless runs
stuck forever - and it also removes the import scratch directory, or the leftover TSVs trip the copy-CSV overwrite guard
(#1596).

`annotation/tasks/annotation_scheduler_task.py:reclaim_stalled_annotation_runs` keys on "occupies a slot but has no live
lease", not only on an expired lease: a worker SIGKILLed before the `finally` leaves a NULL lease that is never `< now`,
and those runs once held slots forever. Disk is a scheduling input too (`annotation/tasks/annotation_scheduler_task.py:has_free_disk_for_annotation`,
#1670): a run's on-disk output is removed by `annotation/signals/annotation_run_cleanup.py` only after its rows commit,
so a stalled import keeps its annotated VCF for the resume lane.

Archived sub-versions (`library/django_utils/data_archive_mixin.py:DataArchiveMixin`) make
`annotation/annotation_version_querysets.py:get_variant_queryset_for_annotation_version` raise `DataArchivedError`;
check `data_archived` before assuming a partition has rows, and restore per `claude/runbooks/restore_partition_archive.md`.
The OntologyVersion is part of the AnnotationVersion, so an ontology re-import invalidates gene annotation on every
build until `gene_annotation` is re-run; `promote_to_active` will tell you, but the variant page would otherwise 500.
