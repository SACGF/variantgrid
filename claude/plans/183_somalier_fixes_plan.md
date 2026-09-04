# Somalier fixes (#183, #162, #196, #393, #432, #1147, SACGF/variantgrid_com#78)

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-04

## Diagnosis in one paragraph

Somalier re-genotypes every sample from the `AD` FORMAT field whenever the VCF header declares `AD`,
and only uses `GT` when the header has no `AD` line. The VCF we hand to `somalier extract`
(`snpdb/variants_to_vcf.py:vcf_export_to_file`) always declares `AD` and writes depths that are
wrong for a large share of real inputs: `.,alt` when the source VCF has `AD` but no `DP` (the AFHCS
files), `ref,-1` / `2×DP,alt` when a depth or AF is stored as the `-1` missing sentinel
(`CohortGenotype.MISSING_NUMBER_VALUE`, which the `is not None` checks let through), and `.` when the
source VCF has no `AD` at all. Reproduced locally on the 166-sample lipo VCF: removing DP the way the
code does for AFHCS takes somalier's het count from 237 to 0 per sample and related pairs from
13,695 to 0. On top of that the `--unknown` decision keys off "all samples came from one VCF", which
is wrong for merged VCFs that contain only variant calls (lipo: 7.3M `./.`, zero `0/0`; without
`--unknown` every pair came out related, median 0.65), the all-samples relate has had no caller since
it was removed from the import pipeline (#393), ancestry crashes on a VCF with zero usable sites and
takes the cohort relate down with it, and `somalier_vcf_id` catches only `CalledProcessError`.

Somalier stage timings measured locally (166 samples): DB export 1.8s for a 3.9M-variant VCF,
`extract` 0.2s, `relate` 0.1s, `ancestry` 34s (it reads the 2,504 1kg `.somalier` files each run).
Ancestry is the cost in #1147.

## Data

No model changes. The data that changes is the VCF we generate, the settings block, and the beat
schedule.

### Somalier extract VCF (generated per VCF, `import_processing/somalier_vcf_extract_<pk>/`)

```
##fileformat=VCFv4.2
##contig=<ID=1,length=248956422,assembly=GRCh38>        # contig names, matching the records and the sites file
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  <slug>_<sample_pk> ...
1       965125  <variant_pk>  G  C  .  .  .  GT  0/1  ./.  1/1
```

`GT` is the only FORMAT field. Somalier then genotypes from `GT` (its own fallback: 20,0 / 10,10 /
0,20 pseudo-depths), so the calls VariantGrid imported and QC'd are the calls somalier relates on.
Verified locally: a GT-only export of HG001 gives somalier exactly the 6,723 hets stored in the DB.

### Settings (`variantgrid/settings/components/default_settings.py`)

```python
SOMALIER = {
    "enabled": False,
    "admin_only": False,
    "vcf_base_dir": ...,
    "report_base_dir": ...,
    "annotation_base_dir": ...,
    "annotation": {...},                 # unchanged
    "ancestry_enabled": True,            # new: the expensive stage, switchable per deployment (#1147)
    "min_genotyped_sites": 100,          # new: a VCF whose best sample has fewer het+hom sites is
                                         #      extracted but ancestry/relate are recorded as SKIPPED
    "all_samples_relate_hour": 2,        # new: nightly all-vs-all relate (#393); None disables
    "relatedness": {                     # unchanged values; the help page reads from here (#183)
        "min_relatedness": 0.1,
        "min_shared_hets": 1000,
        "min_shared_hom_alts": 200,
    },
}
```

### Beat schedule (`variantgrid/celery.py`)

```python
if settings.SOMALIER["enabled"] and settings.SOMALIER["all_samples_relate_hour"] is not None:
    app.conf.beat_schedule["somalier-all-samples-relate"] = {
        "task": "snpdb.tasks.somalier_tasks.somalier_all_samples",
        "schedule": crontab(hour=settings.SOMALIER["all_samples_relate_hour"], minute=0),
    }
```

Route `snpdb.tasks.somalier_tasks.somalier_all_samples` to `SCHEDULING_SINGLE_WORKER` in
`celery_settings.py` so two nightly runs can never overlap.

## §1 Export genotypes somalier will trust (#183)

`snpdb/variants_to_vcf.py:vcf_export_to_file` (somalier is its only caller):

1. Write `FORMAT=GT` and the genotype string only. Keep the per-sample zygosity counting exactly as
   is, since `SomalierSampleExtract` is populated from it.
2. Build the header with `get_vcf_header_from_contigs(..., use_accession=False)` so `##contig` IDs
   match the record `CHROM` values and the `nochr` sites file.
3. Select only the columns the export now uses (`id`, contig name, position, ref, alt, zygosity).
4. Keep the `filters__isnull=True` restriction (somalier only counts PASS sites).

## §2 Decide `--unknown` from the data, not the file count (#183)

`SomalierRelate.is_joint_called_vcf` becomes `has_hom_ref_calls`: true when every sample in
`get_samples()` has a `SomalierSampleExtract` with `ref_count > 0`. A VCF that records hom-ref
genotypes has real `0/0` calls, so an absent site is unknown; a VCF without them (merged
single-sample calls, GIAB benchmark, gVCF-derived variant-only files) means absent is hom-ref, and
`--unknown` is the right flag. `SomalierAllSamplesRelate` always passes `--unknown`, as now.

## §3 Make `somalier_vcf_id` survive anything (#432, variantgrid_com#78, #162)

In `snpdb/tasks/somalier_tasks.py`:

1. Guard against a concurrent run for the same VCF: if a `SomalierVCFExtract` for this VCF is in
   `PROCESSING` and was modified within the last hour, log and return. Otherwise delete it and
   create the new one, as now.
2. Wrap the whole body in `try/except Exception`: store the traceback on the extract
   (`error_exception`), set `ERROR`, log, and return normally. The import step stays green; the
   sample page shows the error (§6).
3. Run each stage independently so an ancestry failure still lets the cohort relate run:
   extract → ancestry (own try) → cohort relate (own try).
4. Skip ancestry and cohort relate with status `SKIPPED` when
   `max(het_count + hom_count)` across the VCF's sample extracts is below
   `SOMALIER["min_genotyped_sites"]`, and skip ancestry when `SOMALIER["ancestry_enabled"]` is
   false. `execute()` already records `PROCESSING/SUCCESS/ERROR`; add the `SKIPPED` writes next to
   the existing `has_genotype` skip.
5. Log a timing line per stage (`somalier <stage> vcf=<pk> took <s>s`) so #1147 can be answered
   from the celery log on prod.

## §4 Bring the all-samples relate back on a schedule (#393, #196)

`somalier_all_samples`:

1. Early exit with `SKIPPED` when no `SomalierVCFExtract` has been modified since the last
   `SUCCESS` run and no sample has been deleted since (compare counts of `.somalier` files against
   `SomalierRelatePairs` sample ids).
2. Build the filtered pairs DataFrame as now, then replace the row-by-row `update_or_create` with:
   delete pairs whose `relate` is not this run, resolve sample ids in one query and drop rows whose
   sample no longer exists, `bulk_create` the rest. On success delete previous
   `SomalierAllSamplesRelate` rows (and their report directories via the existing pre_delete hook).
3. Register the beat entry and routing from the Data section.
4. `somalier_existing_vcfs` keeps calling it directly at the end, as now.

## §5 Help and thresholds (#183, #162)

1. `related_samples_help.html` describes the thresholds in words and the sample page renders the
   live numbers from `settings.SOMALIER["relatedness"]` next to the table heading.
2. Document `--unknown` in the help: relatedness between samples from different VCFs, or from a
   VCF without hom-ref calls, is understated relative to a jointly called VCF (the HSS2008 trio
   measured 0.51 joint vs 0.33 with `--unknown`), so duplicates still show near 1.0 but
   parent/child pairs may sit in the 0.3–0.5 band.
3. Add a wiki page `Install-Somalier.md`: binary + `sites.*.vcf.gz` + `1kg-somalier` +
   `ancestry-labels-1kg.tsv` under `SOMALIER["annotation_base_dir"]`, the sites VCFs imported with
   their file name unchanged (the `deployment_check` command verifies this), the `enabled` flag,
   the `somalier_existing_vcfs --clear` backfill, and a note that local disk matters for ancestry.

## §6 Show what happened on the sample and VCF pages (#162, #196)

`view_sample.html` Ancestry/Relatedness tab and `view_vcf.html`:

1. Show the status and, when present, `error_exception` of the VCF's extract, ancestry run and
   cohort relate, so a failed or skipped stage is visible instead of the tab silently missing.
2. On the sample page list related samples with the patient link state: same patient, different
   patient, or no patient, with a warning badge when relatedness ≥ 0.9 and the patients differ
   (#196 first bullet). A dedicated "duplicate / related samples" page is a follow-up once pairs
   are being produced nightly.

## §7 Deployment

1. A new `snpdb` migration with
   `ManualOperation(task_id=ManualOperation.task_id_manage(["somalier_existing_vcfs", "--clear"]), test=...)`
   where `test` returns true when the deployment has any `SomalierVCFExtract`, so every deployment
   regenerates extracts, ancestry and cohort relates from GT-only VCFs and rebuilds the pairs.
2. `somalier_existing_vcfs --clear` should also delete `SomalierAllSamplesRelate` rows so old
   report directories go with them.

## §8 Tests

Worth keeping:

1. `vcf_export_to_file` writes `GT` only, contig names in header and body, `./.` for unknown
   zygosity, and returns zygosity counts per whitelisted sample (build a small VCF with the
   `snpdb/tests/utils` fixtures; a hom-ref, het, hom-alt and unknown call per sample).
2. `has_hom_ref_calls`: true for a VCF whose extracts have `ref_count > 0`, false when any is 0.
3. `somalier_vcf_id` with the somalier binary replaced by a failing command records `ERROR` on the
   extract and returns without raising; with a fresh `PROCESSING` extract it returns early.
4. `somalier_all_samples` pair loading from a fixture `somalier.pairs.tsv`: threshold filter,
   deleted-sample rows dropped, previous pairs replaced.

## Order of work

§1 and §2 together (they change what somalier sees), §3, §4, §7, then §5 and §6. Re-run
`somalier_existing_vcfs --clear` on vgtest2 first and check the HSS2008 trio and sample 1824
(#183) before rolling to variantgrid.com.
