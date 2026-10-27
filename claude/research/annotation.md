# annotation — research notes

Verified against 7c4408c62 on 2026-09-06

The annotation app takes the Variants snpdb has stored and attaches what the world knows about them: VEP's per-transcript
consequences and the plugin/custom scores that ride along with it, ClinVar's summary per variant, per-gene annotation
(ontology terms, GenCC gene-disease strength, dbNSFP gene scores, gnomAD constraint), Human Protein Atlas tissue data,
and citations. Every one of these is versioned, and `annotation/models/models_annotation_version.py:AnnotationVersion` is the bundle per
build that an analysis or a variant page pins itself to. This document is the story behind `annotation/AGENTS.md`:
how a version comes to exist, how a variant gets its rows, why the tables are partitioned, and the traps.
Model fields, URLs, commands, tasks and signals are in the generated maps
([models](../maps/models.md#annotation), [urls](../maps/urls.md#annotation), [commands](../maps/commands.md),
[tasks](../maps/tasks.md#annotation), [settings](../maps/settings.md)); `claude/domain.md` has the vocabulary.

## Flows

### A new VariantAnnotationVersion

A `annotation/models/models_version.py:VariantAnnotationVersion` (VAV) is a database snapshot of one VEP install for one build:
the VEP code version, the cache version, and a pin for every data file and plugin the command line will reference
(gnomAD, dbNSFP, SpliceAI, COSMIC, MaveDB, Open Targets ...), plus `columns_version`, which selects the set of columns
the pipeline writes. `annotation/vep_annotation.py:get_vep_version` runs the full VEP command over a one-line fake VCF
(cache load included, tens of seconds; `lru_cache`d per process) to read its `##VEP=` header, and
`annotation/vep_annotation.py:vep_dict_to_variant_annotation_version_kwargs` turns that plus
`annotation/vep_config.py:VEPConfig` (the `settings.ANNOTATION[build]` block) into the kwargs;
`annotation/annotation_versions.py:get_or_create_variant_annotation_version_from_current_vep` then `get_or_create`s the
row with `status=NEW`, which is what `manage.py create_new_variant_annotation_version` does per build. Nothing is
re-annotated yet: NEW is a version being built.

The VAV must know which gene set VEP used, because the import lane links every annotation row to a
`TranscriptVersion` by FK. `annotation/gene_annotation_release_matching.py:VepGeneSetVersions` reads the release
token out of the VAV's `refseq` / `genebuild` / `vep_cache` pins (eg `RS_2025_08`), and
`annotation/models/models_version.py:VariantAnnotationVersion.link_gene_annotation_release` looks for an existing
`GeneAnnotationRelease` whose cdot import URL names that release; if there is none, the genes app's
`import_cdot_gene_annotation_release` command downloads the matching cdot per-GFF file and creates it (see
`claude/research/genes.md`). A VAV with no release is a VAV whose every run will fail at pipeline time, which is what
`annotation/models/models_version.py:VariantAnnotationVersion.get_annotation_run_blocker` reports on the runs page. Gene
annotation is then built for that release (`manage.py gene_annotation --new-releases`, below), and only then can
`annotation/models/models_version.py:VariantAnnotationVersion.promote_to_active` run: it refuses while
`get_gene_annotation_promote_blocker` says the AnnotationVersion would be inconsistent, then demotes the prior ACTIVE
VAV for the build to HISTORICAL inside a `select_for_update`. A partial unique constraint (`one_active_vav_per_build`)
makes "the VAV" a database fact rather than a convention.

Saving any sub-version creates its partition tables (`annotation/models/models_version.py:SubVersionPartition.save` via
`library/django_utils/django_partition.py:RelatedModelsPartitionModel`) and calls
`annotation/models/models_annotation_version.py:AnnotationVersion.new_sub_version`, which assembles a new AnnotationVersion from the latest
sub-version of each kind for the build - overwriting the previous AnnotationVersion in place if nothing references it
yet, because a deploy that imports ClinVar, HPA and a new VAV in one go should not leave three half-built versions
behind. `SubVersionPartition.defer_new_sub_version` holds the bump back while a batch of sub-versions is written.
`AnnotationVersion.validate` is the consistency check the rest of the system relies on: every sub-version present, the
GeneAnnotationVersion built for the same GeneAnnotationRelease as the VAV, and the same OntologyVersion.

### Range locks, runs and the two lanes

Annotation is scheduled by variant pk range. `annotation/tasks/annotation_scheduler_task.py:annotation_scheduler`
(beat, and kicked by `upload/tasks/vcf/import_vcf_tasks.py:CheckStartAnnotationTask` after a VCF's variants are
inserted) runs on `scheduling_single_worker` (one process, which is what serialises it with the dispatcher) and, for each
annotated build's latest VAV of the requested status, creates any runs missing from existing locks once
(`annotation/tasks/annotation_scheduler_task.py:_create_missing_runs`), then loops `annotation/tasks/annotation_scheduler_task.py:_handle_variant_annotation_version` until
`annotation/annotation_versions.py:get_annotation_range_lock_and_unannotated_count` finds no more unannotated variants.
That function walks pk blocks above the highest existing lock (never a sort-and-limit over the whole table, which does
not scale) and returns an `annotation/models/models_annotation_run.py:AnnotationRangeLock` of at least
`ANNOTATION_VEP_BATCH_MIN` and at most `ANNOTATION_VEP_BATCH_MAX` variants. "Unannotated" is defined once, in
`annotation/annotation_version_querysets.py:get_variants_qs_for_annotation`: variants in the build's contigs passing
`annotation/models/models.py:VariantAnnotation.get_variant_annotation_q_list` (no reference-alt, no `.`/`*` alt) with no
VariantAnnotation row in this version's partition, restricted to a pipeline type by
`annotation/annotation_pipeline_routing.py:pipeline_type_variant_q`.

Each lock gets one `annotation/models/models_annotation_run.py:AnnotationRun` per enabled pipeline type
(`annotation/tasks/annotation_scheduler_task.py:_handle_range_lock`). The types are the registry in
`annotation/pipelines/__init__.py:PIPELINES`: STANDARD (VEP over short variants), STRUCTURAL_VARIANT (VEP over symbolic
variants from the VAV's `structural_variant_min_size`, and explicit ones whose longer allele reaches its
`explicit_structural_variant_min_size`, #1708), GENE_LEVEL (fusions, computed locally by `annotation/gene_level_annotation.py:annotate_gene_level_run`) and
ANNOTSV, which depends on the SV run having FINISHED and updates the rows it wrote. The scheduler only creates pending
state. `annotation/tasks/annotation_scheduler_task.py:dispatch_annotation_runs` is the single authority that launches:
it sweeps the latest NEW and ACTIVE VAV per build, reclaims stalled runs, fills the count lane
(`annotation/tasks/annotation_scheduler_task.py:count_annotation_runs` pre-counts each run so an empty one finishes
without a dump), merges consecutive pending locks into bigger batches once workers are busy
(`annotation/annotation_versions.py:merge_pending_range_locks`, the inverse of
`annotation/tasks/annotation_scheduler_task.py:subdivide_annotation_run_range`), and leases up to the free capacity
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
`annotation_run_complete_signal`. Status is never set directly: `annotation/models/models_annotation_run.py:AnnotationRun.get_status`
derives it from the timestamp and error fields on every save, and `annotation/models/models_annotation_run.py:AnnotationRun.save_if_owner`
is the write used after any stage during which the lease could have expired. The import lane counts VEP's
`--skipped_variants_file` against the dump so truncated output fails the run (#1701), and SV conservation is scored with
pyBigWig into a sidecar file (#1657).

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
`add_calculated_transcript_columns`. A column added later can be filled from an annotated VCF without re-running the pipeline
(`annotation/backfill_columns.py:import_backfill_vcf`, #1675). `bulk_insert` writes CSVs and `COPY`s them straight into the VAV's partition
tables; under `UNIT_TEST` it falls back to `bulk_create` into the base table, which is why tests never exercise the
partition path.

Reading goes the other way through one hook. Partitions are old-style inheritance children, so a plain ORM join to
`annotation_variantannotation` hits the empty parent; `annotation/annotation_version_querysets.py:get_queryset_for_annotation_version`
attaches `annotation/models/models_annotation_version.py:AnnotationVersion.sql_partition_transformer` to the queryset via
`library/django_utils/django_queryset_sql_transformer.py:get_queryset_with_transformer_hook`, and the compiled SQL has
its table names rewritten to the version's partitions. Django's `FilteredRelation` cannot do this job because it cannot
follow nested relations (`variantannotation__gene__geneannotation`).

### A new column on an existing VariantAnnotationVersion

The default way to add a VEP plugin or data file is a new VAV and a full re-annotation, which takes weeks on a
production database. The alternative keeps the ACTIVE VAV, turns the source on for new variants, and backfills the
variants already annotated. Four pieces have to agree.

**The pin.** Every data file `get_vep_command` passes is pinned on the VAV, either through a dedicated field filled
from the `##VEP=` header or through `annotation/vep_config.py:vep_component_version_kwargs` (basename or parsed
version, from settings alone). A new source means a new VAV field (an entry in `_BASENAME_PIN_FIELDS` /
`_VERSION_PIN_FIELDS`, or header parsing in `vep_dict_to_variant_annotation_version_kwargs`). Setting the file in
`settings.ANNOTATION[build]["vep_config"]` while the ACTIVE VAV still holds NULL in that field is a pin mismatch (see
Traps): the VEP lane halts on that VAV, and `create_new_variant_annotation_version` would `get_or_create` a fresh
VAV. So the field is stamped on the existing NEW/ACTIVE VAVs from live settings at the same time the file is
configured: a `RunPython` that, per VAV, reads `settings.ANNOTATION.get(vav.genome_build_id)`, skips a build with no
config, and sets the new field from what `vep_component_version_kwargs` returns. HISTORICAL VAVs stay NULL - nothing
looks them up by config again. A source that is off in the package defaults
(`variantgrid/settings/components/annotation_settings.py`) leaves the upgrade a no-op; each deployment then turns it on
with a command that stamps the same way from its live settings (`enable_vep_plugin_backfill` below), followed by an
annotation_workers restart, since `get_vep_version` is `lru_cache`d per process.

**The column gates.** The `annotation/vep_columns.py:VEPColumnDef` takes `min_columns_version` equal to the existing
VAVs' `columns_version` (a bump is by definition a new VAV) and `min_vep_version` for the first VEP release whose
VEP_plugins ships the plugin, plus `genome_builds` if the data only exists for some builds. The plugin's settings key
goes in `_VEP_PLUGIN_SETTING_KEYS` so `filter_for(vep_config=...)` drops the columns where the file isn't configured.
`annotation/backfill_columns.py:resolve_backfill_columns` filters by the VAV's own `columns_version` and `vep`, so a
column gated past them cannot be backfilled.

**The backfill flag.** Until the backfill finishes, a NULL in the new column means either "not backfilled" or "the
source has nothing", and for a sparse source those look the same. A `backfilled_<thing>` BooleanField on the VAV
(`annotation/models/models_version.py:VariantAnnotationVersion.backfilled_ptc` and its siblings) separates them:

- `default=True`, because a VAV created later gets the column at insert time.
- The schema migration adds the field, then a `RunPython` sets False only on the VAVs the backfill will actually fill:
  status NEW/ACTIVE, a build the column defs apply to, `columns_version` / `vep` the defs cover. Setting every VAV False
  disables correct data and points operators at work that doesn't exist.
- The command that fills a partition sets the flag True on that VAV when it finishes, and is idempotent so a re-run
  re-asserts it.
- Every reader checks it. An analysis node skips the filter and says so in `get_warnings` on a False VAV
  (`analysis/models/nodes/filters/damage_node.py:DamageNode.get_warnings` for `backfilled_ptc`), since filtering on the
  NULL would drop variants. A display (the variant page) shows nothing, as it does for any missing value: a "not
  backfilled" message there is noise.
- The deploy step is a `ManualOperation` in a migration after the schema one (`manual/__manual_readme.md`), its `test`
  being `VariantAnnotationVersion.objects.filter(<flag>=False).exists()`, with a `note` naming the issue. Use
  `ManualOperation.operation_manage([...])` when one command does the whole job, `operation_other` when a step runs
  outside VariantGrid (external VEP).
- When no VAV is False any more, the field, the readers' fallback branches and the migrations go together.

**The backfill.** `annotation/management/commands/annotation_backfill_columns.py` (#1675, wiki page
"Backfill-Annotation", used for COSMIC counts and Open Targets) dumps already-annotated variants to a VCF carrying
`variant_id`, the operator annotates it outside VariantGrid, and `--import` COPYs the values into a temp table and
UPDATEs the partition. For a VEP plugin: dump bounded by `--max-variant-id` (the last variant annotated before the
plugin was switched on), split by `--min/max-variant-id` to run VEP in parallel, run VEP with the same cache, pick
order and flags as the pipeline so the PICK'd transcript matches, then `--import --csq <column>=<CSQ field>` with
`--dry-run` first. Limits as it stands:

- `--csq` reads the PICK'd CSQ entry into `VariantAnnotation` only; `VariantTranscriptAnnotation` keeps NULL
  (`annotation/backfill_columns.py:_base_table_names_for_column`). A per-transcript column needs the import extended to
  match each CSQ entry to its transcript row, or accepts representative-only history - and a page showing it reads
  `VariantAnnotation`, since the transcript rows of backfilled variants are NULL.
- `--only-missing` dumps rows where every named column is NULL, which for a sparse source never shrinks, so a re-run
  re-dumps the whole version. Bound by variant id instead, and narrow it with `--not-null` (a sibling column the
  source writes) or `--consequence` (an SO term on the representative transcript, all a `--csq` import can fill).
- `--mark-backfilled <flag>` sets the VAV's `backfilled_*` flag once the import finishes. With a split dump, pass it
  on the last import only.

**UTRAnnotator effects (#2499).** Ensembl's UTRAnnotator reported one of a variant's 5' UTR effects, chosen by Perl
hash order; our copy of the plugin reports all of them. `backfilled_utr_effects` is False on the NEW/ACTIVE versions
that ran the old plugin (migration 0009), and `annotation/vep_plugin_backfill.py:backfill_steps` prints one pass per
version for whichever flags it has pending: UTRAnnotator alone over the variants with a consequence
(`--not-null utr_annotator_consequence`) when only this one is, or every source over every variant with no
`--max-variant-id` when the #2426 backfill is pending too, since the variants annotated after that switch-on were
annotated by the old plugin. `--mark-backfilled` takes both flags comma separated. The transcript rows keep the old
plugin's values, as for any `--csq` import.

**UTRAnnotator, AVI and SpliceVault (#2426)** are this recipe applied to three plugins at once, so a deployment runs
one backfill VEP pass and one import:

- UTRAnnotator: five `utr_annotator_*` columns on both annotation tables, GRCh37/38, `min_vep_version=112` (VEP 112
  renamed its uORF stop codon terms).
- SpliceVault: seven `splicevault_*` columns on both tables, GRCh37/38, `min_vep_version=112` (the first VEP_plugins
  release with `SpliceVault.pm`). Its data file is keyed by unversioned ENST and the plugin only returns a row whose
  transcript equals the annotated transcript's `stable_id`, which on a `--refseq` run is a versioned NM - so a
  RefSeq build points `vep_config["splicevault"]` at a copy rewritten to MANE RefSeq ids
  (`annotation/splicevault.py`, `manage.py splicevault_to_refseq`). Rows for non-MANE transcripts are dropped: on
  the start of chr1 that is about 83% of rows (10% of transcripts are MANE). The MANE v1.3 NM versions match the
  VEP RefSeq caches' (215/215 chr21 MANE transcripts in the GRCh38 cache, 214/215 in GRCh37). An Ensembl build uses
  the download as it is.
- AVI (AlphaGenome Variant Impact): `avi_phred` / `avi_raw` on `VariantAnnotation` only (fields declared there
  rather than on `AbstractVariantAnnotation` never reach the transcript rows), GRCh38 SNVs (DeepMind only publishes
  hg38), `min_vep_version=116`. It matches on position and allele, so `--refseq` is fine. `AVI.pm` was added to
  VEP_plugins `release/116` on 2026-08-27, so a VEP 116 plugins checkout made before then lacks it - `vep_data_check`
  reports a configured plugin with no `.pm` in ANNOTATION_VEP_PLUGINS_DIR, and `download_annotation_data` fetches it
  (`annotation/vep_annotation.py:get_vep_plugin_url`).

- AlphaGenome merged splicing (#2433): `alphagenome_splicing` on `VariantAnnotation` only, GRCh38 SNVs,
  `min_vep_version=116`, from a `--custom` VCF (`VEPCustom.ALPHAGENOME_SPLICING`, CSQ field
  `AlphaGenome_splicing_score`). No plugin reads DeepMind's TSV and `--custom` only matches alleles for VCF, so
  `annotation/annotation_data/generate_annotation/alphagenome_splicing_to_vcf.sh` converts it once, with contigs
  stripped of `chr` like our other GRCh38 custom VCFs. The data is licensed for non-commercial use only:
  `ANNOTATION_NON_COMMERCIAL_DATA_ALLOWED` (default True) False clears it in
  `annotation/vep_config.py:allowed_vep_config_data`, which both `VEPConfig` and `vep_component_version_kwargs` go
  through, so a commercial deployment never passes it to VEP, never stamps its pin and never shows the column, even
  with the file configured. The variant page shows a "non-commercial use only" note next to the score.

All are gated `min_columns_version=5` with the `utr_annotator`, `splicevault`, `avi` and `alphagenome_splicing`
basename pins and one flag, `backfilled_utr_avi_splicevault`. `annotation/vep_plugin_backfill.py:SOURCE_PINS` maps each
pin to its plugin or custom source, and only stamps the pins the model it is given has, so an earlier migration's
historical model still works. They are off (None) in the package defaults, so the release changes nothing until a
deployment turns them on. `annotation/vep_plugin_backfill.py:stamp_plugin_pins` stamps the pins on a build's
NEW/ACTIVE versions from live settings and sets the flag False where a pin changed to a file whose columns cover the
version, the version has range locks and it is the build's newest by `annotation_date` (an older one is on its way out,
so isn't worth a backfill); an unchanged pin leaves the flag alone, so it is idempotent.
`annotation/migrations/0004_utr_annotator_avi_splicevault.py` calls it for a deployment that set the files before
migrating, whose ManualOperation (`annotation/migrations/0005_utr_annotator_avi_splicevault_manual_backfill.py`)
carries the steps; `annotation/migrations/0007_alphagenome_splicing.py` does the same for AlphaGenome splicing. Otherwise, per build, under ANNOTATION_VEP_BASE_DIR:

1. Data files. With the paths set in the env settings (step 2), `manage.py download_annotation_data` fetches each
   one, already converted, from variantgrid.com, plus a `.csi` beside each `.tbi` where one is hosted: a 4 kb-bin
   CSI index that htslib prefers over the `.tbi` and that makes VEP's custom/plugin lookups several times faster
   (#2499). A file without a hosted `.csi` is reported and still works on its `.tbi`;
   `annotation/annotation_data/generate_annotation/build_tabix_csi.sh` builds the missing ones in place. The
   sources, for building them from scratch:
   - uORF_5UTR_GRCh37_PUBLIC.txt / uORF_5UTR_GRCh38_PUBLIC.txt from https://github.com/Ensembl/UTRannotator into
     `annotation_data/<build>/` (the plugin computes effects from the transcript sequence and runs with every
     `Evidence` NA without it).
   - SpliceVault_data_hg19.tsv.gz / SpliceVault_data_GRCh38.tsv.gz (+ .tbi, ~880 MB each) from
     https://ftp.ensembl.org/pub/current_variation/SpliceVault/, then for a RefSeq build
     `manage.py splicevault_to_refseq --genome-build <build> --input SpliceVault_data_<hg19|GRCh38>.tsv.gz`, which
     needs MANE loaded and writes `SpliceVault_data_<hg19|GRCh38>.refseq_mane.tsv.gz` + .tbi alongside.
   - GRCh38: "AVI SNV scores" from https://deepmind.google.com/science/alphagenome/downloads (Google sign-in; commercial
     and non-commercial use), `avi_scores_snvs_tabix.zip`, 88.5 GB and (like the splicing zip) stored, so its
     alphagenome_variant_impact_score_snvs.tsv.gz (`#CHROM POS REF ALT raw_score PHRED`, `chr` contigs, which the plugin maps) is about
     as big: extract it and the .tbi into `annotation_data/GRCh38/` from a zip on another filesystem, or stream it
     (`curl <download url> | bsdtar -xf -`) so the zip never lands, as zip + extract need ~177 GB on one disk.
   - `UTRAnnotator.pm` comes from our repo (`annotation/annotation_data/generate_annotation/UTRAnnotator.pm`), not
     Ensembl's: it carries Ensembl/VEP_plugins#869, which reports every 5' UTR effect in a fixed order where the stock
     plugin picks one by hash order and so differs between runs. `download_annotation_data` replaces the installed copy
     whenever it differs from ours (`annotation/vep_annotation.py:is_variantgrid_vep_plugin`). Drop it from
     `_VARIANTGRID_VEP_PLUGINS` once that PR is on the release branch we download from.
   - `AVI.pm`, which `download_annotation_data` fetches from VEP_plugins `release/<ANNOTATION_VEP_VERSION>` into
     ANNOTATION_VEP_PLUGINS_DIR along with any other configured plugin missing there.
   - GRCh38, non-commercial deployments only: "AlphaGenome SNV merged splicing scores",
     combined_splicing_snvs_tabix.zip (20.6 GB, same sign-in), stored rather than compressed: extract
     combined_alphagenome_splicing_snvs.tsv.gz + .tbi, then convert it with
     `annotation/annotation_data/generate_annotation/alphagenome_splicing_to_vcf.sh` to
     annotation_data/GRCh38/alphagenome_splicing_snvs.GRCh38.vcf.gz (20 GB; 42 minutes with 3 jobs on
     vg-test2), after which the TSV can go.
2. Stop annotation_workers, so no VEP run starts on the old settings once the pins are stamped (it would latch
   VEPVersionMismatchError on its version). Env settings: set `ANNOTATION[<build>]["vep_config"]` `utr_annotator`,
   `splicevault` and (GRCh38) `avi` and `alphagenome_splicing` to those paths - the package defaults' comments name them.
3. `manage.py enable_vep_plugin_backfill [--genome-build <build>] [--dry-run]`. It refuses while any file VEP would
   be passed is missing (`annotation/vep_annotation.py:get_missing_vep_files`, which `vep_data_check` prints,
   including a `--plugin` with no `.pm`), and while STANDARD runs on those versions have VEP output not yet imported:
   imported once the inserter expects the new CSQ fields they fail with VEPMissingColumnsError, so let the import lane
   drain or pass `--force` and re-run them. It prints each version's pin changes and backfill steps
   (`annotation/vep_plugin_backfill.py:backfill_steps`, the same text the ManualOperation note holds).
4. Start annotation_workers and restart the other celery workers and gunicorn so they read the settings. If a run
   did latch VEPVersionMismatchError, the re-check on the annotation runs page clears it.
5. The printed backfill:
   1. `annotation_backfill_columns --dump --genome-build <build> --variant-annotation-version <pk> --columns <every
      configured plugin column> --pipeline-type STANDARD --max-variant-id <max> --output vep_plugins_<build>.vcf.gz`
      - every variant the STANDARD pipeline annotated, as AVI scores every SNV (structural and gene-level variants
      stay out). Split by `--min/max-variant-id` to run VEP in parallel.
   2. VEP, the command `annotation/vep_plugin_backfill.py:backfill_vep_command` prints: this build's STANDARD
      pipeline command for the version (so its cache, pick order and PICK'd transcript match), keeping only the
      UTRAnnotator, SpliceVault, AVI and AlphaGenome_splicing `--plugin` / `--custom`, `--fields Feature,PICK,<the CSQ
      fields>`, the output-only flags (`--hgvs --check_existing --af --pubmed --sift --domains --protein --numbers`)
      dropped, `--fork` of the cores less one and `--buffer_size 5000`. The PICK'd entry is the same without those
      flags (all 20,000 of a GRCh38 sample). The time is VEP's own consequence calculation, not the plugins: a sample
      takes as long with none of them as with any one. So send VEP only the variants it can fill (below), not the
      whole dump: on variantgrid.com the whole dump was 45M (GRCh37) and 25M (GRCh38) variants at 100-270 variants/s,
      days of VEP, and the subsets below took hours.
      - UTRAnnotator fills only a 5'UTR PICK and SpliceVault only its own file's positions, so VEP the union of a
        `--consequence 5_prime_UTR_variant` dump and the full dump through `bcftools view -T <SpliceVault chrom/pos,
        'chr' stripped>` (`sort -u` the records under the dump's header): 1.2M of GRCh37's 45M, 1.0M of GRCh38's 25M.
      - AVI and AlphaGenome splicing are position + allele lookups, so GRCh38 skips VEP for them: `bcftools
        annotate` the full dump (main chromosomes only, `bcftools view -t 1,...,22,X,Y`: the dump's header declares
        no contigs and unplaced ones fail) with the AVI TSV (`-c CHROM,POS,REF,ALT,AVI_RAW,AVI_PHRED`, its contigs
        are `chr`-prefixed so `--rename-chrs` there and back) and then, from an indexed file (`-W`), the splicing VCF
        (`--pair-logic exact -c INFO/AlphaGenome_splicing_score:=INFO/score`). One pass reads the whole file, so it
        costs the same for any dump size: about 75 + 110 minutes on one core. It matched VEP on all but 2 of 10,000
        sampled variants, both left-padded MNV-shaped records (`AT>TT`), where bcftools trims to the SNV for AVI and
        `--custom exact` does for the splicing score.
   3. `annotation_backfill_columns --import` with the same `--genome-build` / `--variant-annotation-version`, one
      per file: the VEP subset with the UTRAnnotator/SpliceVault `--columns` and a `--csq <column>=<CSQ field>` each,
      the bcftools file with `--columns avi_phred,avi_raw,alphagenome_splicing --info avi_phred=AVI_PHRED --info
      avi_raw=AVI_RAW --info alphagenome_splicing=AlphaGenome_splicing_score`. `--dry-run` first; the last import
      to finish takes `--mark-backfilled backfilled_utr_avi_splicevault`. The import never touches a variant missing
      from its file, and runs at ~2,400-5,000 records/s. For the VEP file, `--dump-vcf <the VCF sent to VEP>` refuses
      unless the annotated records plus VEP's `--skipped_variants_file` rows add up to its records
      (`annotation/backfill_columns.py:check_annotated_vcf_complete`), as VEP can exit 0 on a short file (#1701).
   4. `rm -f` of the dumps, the annotated VCFs and VEP's warnings and skipped variants files.

The `--csq` import fills `VariantAnnotation` from the PICK'd entry only, so backfilled variants have UTRAnnotator and
SpliceVault values for the representative transcript and NULL transcript rows; SpliceVault values only land when the
representative transcript is the MANE one (VEP's default pick order puts MANE Select first). The variant page reads
them from `VariantAnnotation`. UTRAnnotator reports one `5UTR_consequence` per transcript even when it finds several
effects (the last flagged, in Perl hash order); `5UTR_annotation` keeps every effect.

### ClinVar

ClinVar arrives two ways. The bulk summary is the ClinVar VCF, imported through the upload app's step pipeline
(`annotation/import_task_factories.py` registers the steps): `annotation/tasks/import_clinvar_vcf_task.py:ImportCreateVersionForClinVarVCFTask`
makes a `annotation/models/models_clinvar.py:ClinVarVersion` keyed on the file's sha256 (a re-upload of the same file wipes and
refills the same version's partition), unknown variants are inserted by the normal VCF machinery, then
`annotation/vcf_files/import_clinvar_vcf.py:BulkClinVarInserter` maps INFO fields onto `annotation/models/models_clinvar.py:ClinVar`
rows, including the oncogenic (ONC*) and somatic clinical impact (SCI*) fields (#1791). Only staff may upload one
(`annotation/vcf_files/import_clinvar_vcf.py:check_can_import_clinvar`). The per-submission detail is fetched lazily
from ClinVar's XML: `annotation/clinvar_fetch_request.py:ClinVarFetchRequest.fetch` caches a
`annotation/models/models_clinvar.py:ClinVarRecordCollection` per variation id, refreshing when older than
`CLINVAR_RECORD_CACHE_DAYS` or when `annotation/clinvar_xml_parser_via_vcv.py:ClinVarXmlParserViaVCV` bumps its
`PARSER_VERSION`, holding a `select_for_update` so two page loads do not both hit the API.

### Gene annotation

`annotation/models/models_gene_annotation.py:GeneAnnotation` is a per-gene, per-`annotation/models/models_gene_annotation.py:GeneAnnotationVersion` row
of pipe-joined OMIM/HPO/MONDO terms, GenCC strength flags, a dbNSFP gene row and gnomAD LoF o/e. It is built by
`annotation/management/commands/gene_annotation.py:Command`, which must be re-run whenever the GeneAnnotationRelease or
the OntologyVersion changes because the version is unique on both: `_populate_gene_annotation_version` matches every
symbol into the release with `genes/gene_matching.py:ReleaseGeneMatcher`, refreshes PanelApp Australia in one crawl,
then walks every HGNC term through an `ontology/ontology_traversal.py:get_ontology_traverser` (in-memory graph optional)
and groups the results by gene id before a bulk write. A version that matches nothing raises rather than leaving an
empty version live, since an empty GeneAnnotation partition makes every OMIM/HPO filter silently match nothing.
`--missing` and `--new-releases` are the deploy-time modes; `--force` rebuilds.

## Why it is shaped this way

Three "versions" are kept apart on purpose (`annotation/AGENTS.md` says how to tell them apart). The VAV is what VEP
ran with, and its partition is the unit of retention: old analyses keep reading their old partition after an upgrade,
and dropping a HISTORICAL VAV drops one table rather than deleting millions of rows. `columns_version` exists so that a
VEP upgrade that changes nothing in the schema is a new VAV without a migration, and a schema change (new plugin
columns) is a new `VEPColumnDef` gated by `min_columns_version`; `annotation/vep_columns.py:VEP_COLUMNS` is code, not a
table, so the mapping is versioned with the code that reads it and can be tested and diffed. `columns_version` 2 adds
dbNSFP rankscores, 3 gnomAD v4 and the `_xy` fields, 4 raw dbNSFP 5.3 scores with per-transcript resolution
(`annotation/refseq_ensembl_resolver.py`), denovo-db and the VEP 116 plugins (#1652). `annotation/models/models_annotation_run.py:AnnotationPipelineVersion`
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

External annotation (#1568, `annotation/external_annotation.py`) is for deployments where VEP is too slow on the web
host: runs on a NEW VAV can be dumped as self-describing gzipped VCFs with a sidecar of the VAV identity and
the exact VEP command templated for another machine, and imported back into the normal import lane. The range endpoints
are verified against local Variant pks on the way back (`annotation/external_annotation.py:verify_annotated_vcf_variant_ids`)
because the dump's `variant_id` INFO is the only link between the external file and the database.

## Traps

The pins on a VAV must match the VEP that will run and the VCF that comes back: `vep_check_command_line_version_match`
and `annotation/vep_annotation.py:vep_check_annotated_file_version_match` both raise `VEPVersionMismatchError`, so
changing a `settings.ANNOTATION` data path without creating a new VAV halts every run. `gencode_subset` and `distance`
are settings snapshots, not VEP header facts, and are excluded from the comparison
(`annotation/vep_annotation.py:_vep_check_version_match`). The first mismatch latches the VAV
(`annotation/models/models_version.py:VariantAnnotationVersion.latch_vep_version_mismatch`, #1679), which stops the scheduler and the
dispatcher's VEP lane for it; only the runs page's re-check button
(`annotation/tasks/annotate_variants.py:recheck_vep_version`) clears it.

`annotation/models/models_annotation_run.py:AnnotationRun.get_for_variant` answers "is this variant in progress" by range-lock bounds
and pipeline type, so a variant inside a lock with no row reads as pending, not missing; an unannotated count must apply
`get_variant_annotation_q_list` or it will count reference variants that never get rows. Range locks and ClinVar rows PROTECT
their variants (`annotation/models/models_annotation_run.py:AnnotationRangeLock.release_variant` is the way to move a lock's endpoint
off a variant you need to delete).

Merging is only safe while every run on a lock is still CREATED and un-leased
(`annotation/annotation_versions.py:_range_lock_is_dispatchable`); an empty-finished run is reopened when its lock grows
because it owns no rows (`annotation/models/models_annotation_run.py:AnnotationRun.reopen_to_created`). A retry must reuse the run in
place (`annotation/models/models_annotation_run.py:AnnotationRun.reset_for_retry`, #1654) - deleting and recreating it leaves a rangeless
run stuck forever - and it also removes the import scratch directory, or the leftover TSVs trip the copy-CSV overwrite guard
(#1596).

`annotation/tasks/annotation_scheduler_task.py:reclaim_stalled_annotation_runs` keys on "occupies a slot but has no live
lease", not only on an expired lease: a worker SIGKILLed before the `finally` leaves a NULL lease that is never `< now`,
so a check on expiry alone would leave that run holding its slot forever. Disk is a scheduling input too (`annotation/tasks/annotation_scheduler_task.py:has_free_disk_for_annotation`,
#1670): a run's on-disk output is removed by `annotation/signals/annotation_run_cleanup.py` only after its rows commit,
so a stalled import keeps its annotated VCF for the resume lane.

Archived sub-versions (`library/django_utils/data_archive_mixin.py:DataArchiveMixin`) make
`annotation/annotation_version_querysets.py:get_variant_queryset_for_annotation_version` raise `DataArchivedError`;
check `data_archived` before assuming a partition has rows, and restore per `claude/runbooks/restore_partition_archive.md`.
The OntologyVersion is part of the AnnotationVersion, so an ontology re-import invalidates gene annotation on every
build until `gene_annotation` is re-run; `promote_to_active` will tell you, but the variant page would otherwise 500.
