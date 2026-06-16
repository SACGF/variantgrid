# paper/scripts — read-only usage-mining for the VariantGrid paper

**Claude does not run these. You do, after review, on a non-clinical-safe path (read replica or
restored backup preferred).** See `../data_collection_plan.md` §0 for the full safety contract.

## What's here
- `paper_data_mining.py` — a Django management command that emits **aggregate-only** CSVs
  (counts/rates/distributions), with small-cell suppression and no identifiers.

## Run
```bash
# 1. copy into any app's management/commands dir, e.g.:
cp paper/scripts/paper_data_mining.py snpdb/management/commands/paper_data_mining.py
# 2. run (read-only; aggregate CSVs only)
python manage.py paper_data_mining --output-dir /tmp/vg_paper_stats --min-cell 5
#    options: --skip-heavy   (skip slow collectors)
#             --only ingestion_over_time,node_type_usage
```

## Guarantees built in
- Read-only (counts/aggregations only — no writes/migrations).
- No PII: no patient IDs, sample names, usernames, free-text. Users = distinct *counts* only.
- Cells with count < `--min-cell` are suppressed (`<5`).
- `_run_manifest.csv` records what ran, what was emitted, and any failures (auditable provenance).

## Before you trust the numbers
- Run after migrating `vg3_sapath_prod` data up to VG4 (scripts target current models).
- Each `# VERIFY` comment marks a field/relation that may differ across schema versions — confirm
  in-context (a couple relate to `VCF.date`, the diff per-column field, and the recurrence join).
- Sanity-check 2-3 totals against known VG3 values to confirm the migration preserved them.
- Collectors fail independently — one error won't stop the others; check `_run_manifest.csv`.

## Outputs (one CSV each) → see mapping in `../data_collection_plan.md` §2
ingestion_by_month · genotype_variant_ratio · annotation_versions · annotation_version_diffs
(+ _diff_columns) · classifications_by_month (light) · significance_changes_by_month ·
analyses_by_month · node_type_usage · template_runs_by_month · variant_tags_by_type / _by_month ·
classified_variant_recurrence.
