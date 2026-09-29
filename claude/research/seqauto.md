# seqauto — research notes

Verified against a96540a68 on 2026-09-25

seqauto is the lab's sequencing ledger: which flowcells were run (`seqauto/models/models_seqauto.py:SequencingRun`),
what the sample sheet said was on them (`SampleSheet`, one `SequencingSample` per row), the files the pipeline made
from each row (FastQs, BAM, per-sample and joint-called VCFs), the QC measured along the way, and enrichment-kit gold
coverage. It does not find those files itself any more - the sequencing pipeline posts records into `/seqauto/api/`
and VCF / coverage files through the upload API, and seqauto's job is to join the two by path so an imported `Sample`
knows which sheet row, run and kit it came from. It only switches on where `SEQAUTO_ENABLED` (SA Pathology); the
DRAGEN TSO 500 pair records (`LibraryQC`, `DragenTSO500CombinedVariantOutput`) live here because they are keyed on a
run. Fields, URLs, commands, tasks and signals are in the generated maps ([models](../maps/models.md#seqauto),
[urls](../maps/urls.md#seqauto), [tasks](../maps/tasks.md#seqauto), [commands](../maps/commands.md),
[signals](../maps/signals.md)); `seqauto/AGENTS.md` has the short rules.

## Flows

### The pipeline posts a run

The client (`seqauto/management/commands/seqauto_api_client2.py` is the reference, `seqauto_api_client.py` the
hardcoded test-data version) posts in dependency order. `seqauto/serializers/sequencing_serializers.py:SequencingRunSerializer`
`get_or_create`s by name - the run's name is its primary key, so a re-post is a no-op and never updates fields.
`SampleSheetSerializer.create` `update_or_create`s the sheet by (run, hash), upserts its rows by `sample_id`, stores
the extra sheet columns as `SequencingSampleData` key/values, carries extraction links forward from the previous
current sheet (`carry_extractions_to_new_sample_sheet`), then calls
`seqauto/models/models_seqauto.py:SampleSheet.set_as_current_sample_sheet`: whatever sheet was posted last is current.
`SequencingRunCurrentSampleSheet` is the one-to-one pointer that lets queries walk run -> current sheet, and
`SequencingSample.get_current` is "rows on a current sheet".

Per-sample files come in one bulk call, `seqauto/views_rest.py:SequencingFilesBulkCreateView` ->
`SequencingFilesBulkCreateSerializer`, which resolves each record's `sample_name` on the sheet and creates
UnalignedReads (optional since the FastQ-less pipelines, SACGF/variantgrid#357), BamFile and `SingleSampleVCF`.
Joint calls post through `JointCalledVCFSerializer`; a joint call spanning runs (a trio sequenced on different
flowcells) sends explicit `sequencing_samples`, and an empty M2M means "the whole sheet"
(`seqauto/models/models_seqauto.py:JointCalledVCF.get_sequencing_samples`). QC hangs off a `QC` row found by
(sequencing sample, BAM path, VCF path) in `seqauto/serializers/seqauto_qc_serializers.py:QCSerializer.get_object`,
with `QCExecSummary`, `QCGeneList` and `QCGeneCoverage` bulk endpoints beside it. Writes need a superuser or
`SEQAUTO_API_WRITE_GROUP` (`seqauto/views_rest.py:SeqAutoWritePermission`); reads are open.

### The VCF arrives and is joined by path

The VCF file itself is uploaded separately with `path` set to the same string the pipeline gave the record.
`upload/vcf/vcf_import.py:create_backend_vcf_links` looks the path up as a `JointCalledVCF`, then a
`SingleSampleVCF`, and makes a `BackendVCF`; a path matching neither fails the import, which is why
`seqauto/serializers/sequencing_serializers.py:validate_unique_vcf_path` refuses the same path on two records.
`upload/vcf/vcf_import.py:link_samples_and_vcfs_to_sequencing` then copies what seqauto knows onto the VCF and its
samples: zygosity-count flag and `variants_type` from the kit, fake-data flag from the run, the BAM path as a
`SampleFilePath`, the sequencing sample's extraction, a `SampleFromSequencingSample` link, a `VCFFromSequencingRun`
row for the run page (a re-analysis from the same caller replaces the arm's earlier one), pending QC gene lists as
`SampleGeneList`s, and any TSO 500 analysis waiting on that arm (`DragenTSO500CombinedVariantOutput.link_arm_sample`).

Matching VCF samples to sheet rows is `seqauto/models/models_seqauto.py:get_samples_by_sequencing_sample`: upper-case,
`-` to `_`, and accept either the sheet name or bcl2fastq's `<name>_S<number>`; a single-sample VCF whose one sample is
named something else (SpliceGirl's `SAMPLE`) falls back to the filename, longest match wins. No match raises, and the
error lists the names it would have accepted.

Gene coverage works the same way: `upload/tasks/import_gene_coverage_task.py` finds the `QCGeneCoverage` whose path
the upload carries and takes the genome build and kit from it. `QCGeneCoverage`'s `pre_delete` handler deletes the
`GeneCoverageCollection` itself rather than let the cascade run, and swallows the failure when a gold reference
protects it (SACGF/variantgrid_sapath#395).

### A sheet is re-sent

A second sheet for a run becomes current, and `seqauto/sequencing_files/sample_sheet.py:current_sample_sheet_changed`
compares it with the old one (`sample_sheet_meaningfully_changed`: same sample names and barcodes). If nothing
meaningful changed it re-points `SampleFromSequencingSample` links and cross-run joint-call members at the new rows
and moves FastQs / UnalignedReads / BAMs across by sample name. If it did change, nothing moves: the run page shows
"data out of date" (`SequencingRun.is_data_out_of_date_from_current_sample_sheet`) and an admin presses the button
that runs `assign_old_sample_sheet_data_to_current_sample_sheet`, which relinks everything and re-runs
`link_samples_and_vcfs_to_sequencing(replace_existing=True)`. Only the automatic path sends
`sequencing_run_current_sample_sheet_changed_signal`. seqauto itself receives none of its signals; SA Pathology
connects them in `variantgrid_sapath/sapath/apps.py` - the sheet-created one sets kits from the sheet's Panel
columns, the sheet-changed one relinks samples to patients.

### DRAGEN TSO 500 pairs

A TSO 500 pair is a DNA and an RNA library from one specimen. Two uploads describe it, in either order and either
before or after the run is posted: MetricsOutput writes one `LibraryQC` per (run name, pair, category)
(`upload/tso500/dragen_metrics_output_records.py:write_library_qc`) and CombinedVariantOutput writes one
`DragenTSO500CombinedVariantOutput` per (run name, pair) with TMB / MSI / GIS
(`upload/tso500/dragen_combined_variant_output_records.py:write_combined_variant_output`). Both store the run name as
text beside a nullable run FK, because a nullable FK cannot be part of a unique key. Both claim a Specimen by the
accession inside the pair ID without creating one (`seqauto/models/models_seqauto.py:SpecimenClaimMixin`), and each
arm finds its sheet row through the `Pair_ID` / `Sample_Type` sheet columns
(`seqauto/models/models_seqauto.py:sequencing_sample_for_pair`). Whatever cannot be linked yet is parked and swept by
`patients/tasks/extraction_matching_tasks.py:reconcile_pending_extractions`, `link_library_qc_to_sequencing_samples`
and `link_combined_variant_outputs`. `seqauto/views.py:view_tso500_pair` is the pair's page; the lab's MSI / TMB
calls are properties over the `TSO500_*_CALL_BANDS` settings (`DragenTSO500CombinedVariantOutput.msi_call`).

### Gold coverage and stats

`manage.py set_gold_standard_runs` flags runs as gold for a kit and runs
`seqauto/tasks/gold_summary_tasks.py:calculate_gold_summary` in-process (`.apply()`), which snapshots the collections
used (`GoldGeneCoverageCollection`) and bulk-creates per-gene `GoldCoverageSummary` stats. The gold-coverage REST views
and the gene grid's coverage column read these. `seqauto/models/models_seqauto.py:get_20x_gene_coverage` counts samples at 100%
20x for a gene incrementally, caching (count, collections, max pk) for 30 days; classification autopopulate and
`variantopedia/views.py` call it. `seqauto/seqauto_stats.py` feeds the runs-per-month graphs, and
`seqauto/signals/seqauto_integration_status.py` reports the last run / sheet / QC on the integration status page.

## Why it is shaped this way

seqauto began (around 2016, from the TAU project) as a filesystem scanner: a celery beat task ran shell `find` scripts over
the sequencing directories, created a record per file keyed on its path, and wrote PBS job scripts for whatever was
missing. That is why every file model is a `SeqAutoRecord` with `path`, `file_last_modified` and `hash`, why `QC`
is a SeqAutoRecord with a path it does not need (its docstring says so), and why the `get_params` chains and
`SEQAUTO_*_PATTERN` settings exist: the path of the next file was derived by string substitution down
run -> sheet -> sample -> BAM -> VCF. The API (issue #76, 2024) kept those models so the scanner and API could run
side by side, then #1643 deleted the scanner. What survived is the path as the join key between a record the pipeline
posts and a file it uploads - the two arrive through different endpoints, often from different processes, and the
path is the one string both know.

The run name as primary key comes from Illumina run directories being globally unique; it makes the API idempotent
without lookups, at the cost that `SequencingRunSerializer` cannot update a run. Sheets are versioned rather than
edited because labs do re-demultiplex with corrected sheets, and the "meaningfully changed" test is what decides
whether existing BAMs and samples can be trusted to still mean the same library.

`SequencingSample` carries `ExtractionMatchMixin` (#1704) because the extraction is a property of the library, not of
whichever sheet described it - hence the carry-forward on re-send - and `link_samples_and_vcfs_to_sequencing` pushes
it down to every Sample so one link call covers all of an arm's VCFs.

## History

- 2016-2020: TAU-derived scanner; `SeqAutoRun`, `SeqAutoMessage`, `data_state`, job scripts. Public repo from 2020.
- 2024 (#76): REST API and serializers; the API client command.
- 2025: sequencer get-or-create over the API (SACGF/variantgrid_sapath#364); single-sample VCFs matched by filename
  (SACGF/variantgrid_sapath#362); QC gene list duplicates between scan and API (#361, sapath#382); QCGeneCoverage
  cascade (sapath#395).
- 2026-05: `VCFFile` renamed `SingleSampleVCF`, and multi-sample trio VCFs (#1555, #343).
- 2026-07 (#1643, 64dd39081, c1ad599da): the scanner, job scripts, `SeqAutoRun`, `SeqAutoMessage`, `data_state`,
  the InterOp / FastQC / flagstat parsers and their commands deleted; `seqauto/sequencing_files/sample_sheet.py` is
  what was kept of the scanner's create_resource_models module.
- 2026-08: cross-run joint-called VCFs (SACGF/variantgrid_sapath#415); FastQ-less file posts; Specimen -> Extraction
  -> Sample and the sequencing-sample extraction link endpoint (#1704, #1707).
- 2026-09: `LibraryQC`, `DragenTSO500CombinedVariantOutput` and the pair page (#1904); API writes
  gated on `SEQAUTO_API_WRITE_GROUP` (variantgrid_private#3875).

## Traps

- Paths are compared as exact strings. The upload's `path` must be byte-identical to the record's, and one path can
  belong to only one `SingleSampleVCF` / `JointCalledVCF` (`validate_unique_vcf_path`).
- `SequencingRunSerializer.create` is `get_or_create`: re-posting a run with a new kit, experiment or flags changes
  nothing. Edit it on the run page (staff) or admin.
- `seqauto/serializers/sequencing_serializers.py:SampleSheetSerializer.update` is broken: it pops
  `sequencing_samples` (the field is `sequencingsample_set`) and calls `instance.sequencing_samples`, so a PUT/PATCH to
  the sample sheet endpoint raises. Were it fixed as written it would delete the rows and cascade their BAMs and VCF
  records. The pipeline only ever POSTs a new sheet.
- `sequencing_run_created_signal` is never sent since #1643 (the scanner sent it), so
  `sapath_sequencing_run_created_handler` is dead. The sheet-created handler falls back to the kit from the run path,
  so kits still get set when a sheet arrives - but a run posted with no sheet gets no kit.
  `backend_vcf_import_success_signal` is sent with no receiver.
- `seqauto/views.py:assign_data_to_current_sample_sheet` and `reload_experiment_name` check no permission. The first
  relinks a run's data and samples (`replace_existing=True`); its button is on the superuser-only Admin tab, but any
  logged-in user can POST it (compare `delete_sequencing_run`, superuser, and the run form, staff). The second's button
  is on the Experiment tab for everyone, and it reads `RunParameters.xml` from the run's path on the web host, which
  only works where the sequencing filesystem is mounted. The Admin tab's text still says deleted data "will
  re-generate next disk scan" - there is no scan; the pipeline has to re-post.
- `seqauto/models/models_seqauto.py:get_20x_gene_coverage` increments with `gene_coverage_collection__pk__gte` the
  cached max pk, so the collection at the old max is counted again on every increment; and the transcript query is not
  restricted to current-sheet collections the way the collection count is. Numbers drift up over time until the cache
  is invalidated.
- `sequencing_run.save()  # Re-validate ready` in `seqauto/sequencing_files/sample_sheet.py` is a leftover - the
  `ready` field is long gone. Likewise most `SEQAUTO_*` path settings are scan-era and read by nothing; only the QC,
  gene-coverage and GOI patterns are still used, to derive a default path when the API omits one
  (`QC.get_path_from_vcf`, `QCGeneList.get_path_from_qc`). `seqauto/scripts/tau/` is the scanner's shell scripts, unused.
- `SampleFromSequencingSample.sample` is one-to-one but `sequencing_sample` is not: an arm has a Sample per caller
  VCF and per re-import. `replace_existing` only decides whether an already-linked Sample is re-pointed at a new sheet
  row; each newly imported Sample always gets its own link.
- The sheet-created signal fires inside `set_as_current_sample_sheet` before the current-sheet pointer moves, so a
  receiver must use the `sample_sheet` it is given, not `sequencing_run.get_current_sample_sheet()`.
- `get_samples_by_sequencing_sample` normalises only case and `-`/`_`; a VCF sample named anything else raises
  "Couldn't link VCF samples to sequencing samples" and the whole upload fails.
