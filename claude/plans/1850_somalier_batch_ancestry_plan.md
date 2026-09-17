# #1850 — somalier_existing_vcfs: batch ancestry across VCFs, fan extracts out to celery

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-09
Status: draft

[#1850](https://github.com/SACGF/variantgrid/issues/1850). Grew out of the #1147 measurements, which are
the numbers below. #1842 (skip ancestry/relate below a minimum number of genotyped sites) is already in
the code as `settings.SOMALIER["min_genotyped_sites"]` and this plan keeps that behaviour per VCF.

## The problem

`snpdb/management/commands/somalier_existing_vcfs.py` walks every VCF without a `SomalierVCFExtract` and
calls `snpdb/tasks/somalier_tasks.py:somalier_vcf_id` for each, in-process and serially, then runs
`somalier_all_samples(force=True)`. Per VCF (#1147, 166-sample VCF, somalier data on local disk):

| stage | seconds | notes |
|---|---|---|
| DB export (`vcf_export_to_file`) | ~2 | genotypes at the ~17k somalier sites |
| `somalier extract` | 0.2 | |
| `somalier ancestry` | ~34 | fits PCA + a small NN on the 2,504 1kg `.somalier` files - every call, regardless of how many query samples |
| cohort `relate` | 0.1 | |

Ten thousand VCFs is ~100 hours, ~94 of them ancestry, on one core. Nothing is parallel and nothing is
shared between VCFs.

`somalier ancestry --labels L 1kg/*.somalier ++ q1.somalier q2.somalier …` takes any number of query
files after `++`, and somalier expands globs itself (`SomalierAllSamplesRelate.get_sample_somalier_filenames`
already relies on that), so the background fit can be paid once per batch instead of once per VCF.

## Data

Ancestry currently runs against one `SomalierVCFExtract`: `SomalierAncestryRun` is one-to-one with it and
owns the `uuid` that names the report dir under `SOMALIER["report_base_dir"]/ancestry/`. A batch run
produces one report (HTML + TSV) covering every sample in the batch, so the report and its uuid move to a
batch, and the per-VCF run keeps its own status (a VCF can be skipped for too few sites while its
neighbours run) and points at the batch whose report it is in.

```python
class SomalierAncestryBatch(AbstractSomalierModel):
    """ One 'somalier ancestry' call. status/error_exception come from AbstractSomalierModel """
    uuid = models.UUIDField(default=uuid.uuid4, editable=False)  # names the report dir


class SomalierAncestryRun(AbstractSomalierModel):
    vcf_extract = models.OneToOneField(SomalierVCFExtract, on_delete=CASCADE)
    batch = models.ForeignKey(SomalierAncestryBatch, null=True, on_delete=CASCADE)  # null = SKIPPED
```

`SomalierAncestry` is unchanged (FK to run, one-to-one with `SomalierSampleExtract`).

Migration: `AddField batch`, create one `SomalierAncestryBatch` per existing run copying its `uuid` and
status, point the run at it, then `RemoveField uuid` from the run. Report dirs on disk stay valid because
the uuid is preserved. Rows with status SKIPPED / ERROR-before-execute keep `batch = None`.

Report-dir ownership moves with the uuid: `somalier_ancestry_run_pre_delete_handler` becomes a handler on
`SomalierAncestryBatch`. Deleting a run (cascade from `--clear`) leaves the batch report alone; deleting
the batch cascades to its runs and removes the dir.

`SomalierAncestryRun.url` resolves through `batch`, so the Ancestry tab on the VCF page
(`snpdb/templates/snpdb/data/view_vcf_cohort.html`) and `get_stages` keep working with no template
change. The report for a batched VCF shows every sample in its batch against the 1kg background rather
than just its own; that is acceptable - the plot is 2,504 background points either way and the per-sample
prediction is what the sample page shows.

## Settings

```python
SOMALIER = {
    …
    "ancestry_batch_samples": 1000,   # max query samples per 'somalier ancestry' call (backfill only)
}
```

Batched by sample count, not VCF count, so a run of big joint-called VCFs doesn't build a batch of a
hundred thousand samples. Bounded so one bad `.somalier` file fails a batch of a thousand, not the lot.

## Tasks (`snpdb/tasks/somalier_tasks.py`)

Split `somalier_vcf_id` into stages that can be run for one VCF or many:

- `somalier_vcf_extract(vcf_id)` - celery task, queue `db_workers`. Today's extract + zygosity counts,
  then the cohort relate (0.1s, per VCF anyway). Records status on the extract / relate rows exactly as
  now; never raises.
- `somalier_ancestry_batch(vcf_extract_ids)` - celery task, queue `db_workers`. For each extract: create
  the `SomalierAncestryRun`, and set it SKIPPED with the reason if `ancestry_enabled` is off or
  `_max_genotyped_sites` is below `min_genotyped_sites`. If any runs are left: create the batch, run one
  `somalier ancestry` over their `.somalier` files, then `_load_ancestry_tsv(batch, runs)` reads the
  TSV once and creates the `SomalierAncestry` rows per sample extract. Batch status → every member run.
  Query files are passed as `<batch temp dir>/*.somalier`, a directory of symlinks under
  `get_import_processing_dir`, so a batch of a thousand samples never nears `ARG_MAX`; the symlink dir
  goes with the other temp files under `IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS`.
- `somalier_vcf_id(vcf_id)` - unchanged signature and behaviour for the import pipeline
  (`upload/tasks/vcf/genotype_vcf_tasks.py:SomalierVCFTask`): `somalier_vcf_extract` then
  `somalier_ancestry_batch([extract.pk])`. A batch of one, same code path.

`_somalier_ancestry` takes the batch and the list of `.somalier` paths; the per-VCF TSV loop it has
today moves into `_load_ancestry_tsv`, keyed by `AbstractSomalierModel.sample_name` as now.

## The command

```
somalier_existing_vcfs [--genome-build B] [--clear] [--no-wait]
```

1. `--clear` as today (also deletes `SomalierAncestryBatch`, which takes the report dirs).
2. Dispatch `somalier_vcf_extract.delay(pk)` for every VCF with no extract (per build when asked), in pk
   order, logging the count. Existing `_create_vcf_extract` skip-if-still-processing logic stays, so a
   re-run of the command only picks up what is missing.
3. Wait until no VCF is without an extract or has one in PROCESSING younger than
   `STALE_PROCESSING_AGE`; poll the DB every 30s and log progress (done / total). Polling rather than a
   chord because the command is then restartable at any point and needs no result-backend state: kill
   it, run it again, it continues.
4. Group extracts with SUCCESS status and no `SomalierAncestryRun` into batches of
   `ancestry_batch_samples` samples (walking extracts in pk order, summing `somaliersampleextract`
   counts) and dispatch `somalier_ancestry_batch.delay(ids)` for each; wait the same way.
5. `somalier_all_samples.delay(force=True)`.

`--no-wait` dispatches step 2 and exits with a line saying to re-run without it once the queue drains -
for a deploy where the operator would rather not hold a terminal for an hour.

Expected on 10k VCFs / 30k samples with four `db_workers` processes: extracts ~1.5 h (DB-bound, 2s each),
ancestry 30 batches × ~35s ≈ 20 min, all-samples relate ~45 min (somalier: ~1 min for 4,500, O(n²)).
Under 3 hours from ~100.

## Not in scope

- `vcf_export_to_file` cost (2s per VCF). It is the cohort-genotype join at the sites; once extracts
  are parallel it is no longer the wall-clock bottleneck. Measure on the real box before touching it.
- The nightly `somalier_all_samples`: already one call, and somalier drops pairs with relatedness
  ≤ 0.05 by default so `_load_somalier_pairs` stays small.

## Tests (`snpdb/tests/test_somalier.py`)

Keep the ones that cover our logic:

- `_load_ancestry_tsv` splits one TSV with samples from three VCFs into the right `SomalierAncestry`
  rows per run, and a sample in the TSV not in any run is ignored.
- `somalier_ancestry_batch`: a VCF below `min_genotyped_sites` gets a SKIPPED run with `batch=None`
  while its neighbour runs; when every VCF is skipped no batch row and no somalier call.
- A failing `somalier ancestry` marks the batch and every member run ERROR, none raised
  (the existing `test_failure_is_recorded_not_raised` pattern).
- Batching by sample count: extracts of 600, 600, 300 samples with a limit of 1000 give [600], [600, 300].
- Existing `SomalierVCFTaskTest` still passes against the unchanged `somalier_vcf_id`.

The migration's uuid copy is checked by hand on this box (`vg-test2` has 21 extracts): report URLs on
VCF pages resolve before and after.

## Definition of done

- `scripts/vg tests --explain` names `snpdb.tests.test_somalier` and it passes.
- `somalier_tasks.py` docstring lists the three entry points and which is the import path.
- One line in `snpdb/CLAUDE.md` next to the existing somalier note: ancestry is batched, the report a VCF
  links to is its batch's.
- `scripts/vg map` refreshed (new model, task, setting); `scripts/vg docs check` passes.
