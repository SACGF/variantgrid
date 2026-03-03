# VG3 SA Path Production Branch — Commit Summary

**Branch:** `vg3_sapath_prod`
**Diverged from master at:** 2021-06-22 (commit `62035452d` — #170 Haem Mutect auto-load)
**Unique commits (not in master):** 518
**Date range of unique commits:** 2021-06-25 → 2026-01-28

This branch represents the SA Pathology deployment of VG3, receiving cherry-picked fixes and limited backports from master.

---

## Key Features Unique to / Cherry-Picked into VG3

### Sequencing Automation (SeqAuto)
- SAP #360 — Show SeqAuto scan run on grid + scheduled task fix
- SAP #362 — Allow matching on VCF for single-sample VCF
- SAP #364 — NovaSeq X support: get-or-create sequencer
- SAP #365 — Be able to create sequencer via API
- SAP #373 — Auto-run analyses for onco (hookup + signal fix)
- SAP #382 — Fix: Sample GOI duplicated from API/Sequencing scan
- SAP #385 — Onco auto-analysis fixes
- SAP #395 — Better handle reload of gene coverage (stop cascade delete)
- SAP #334 — Prefer specific FASTQ style, only take 1
- SAP #338 — SeqAuto single-sample VCF support
- SAP #353 — VCF link to backend combo/single VCF path - duplicates

### Analysis / Candidate Search
- VG #1428 — Analysis templates: view how samples would match auto runs
- VG #1398 — Output node CSV download (with styling)
- VG #385 — Make auto analysis jobs run via Celery (don't take down whole import)
- VG #1396 — Shrink phenotype terms display

### Infrastructure / Deployment
- SAP #383 — Move VG deployment to new server
- SAP #377 — Limit maximum size of VEP annotation
- SAP #384 — Fix: VariantTag liftover failing
- SAP #357 — Fix: ValueError: None is not a valid ClinVarReviewStatus
- VG #1414 — Fix: Uncaught SyntaxError (Octal escape) in template strings
- VG #1378 — Site messages now allow Markdown
- Seqauto scan disabled by default (still allow manual scans)
- Fix crash bug retrieving ClassificationModificationSync record #1415
- Attempt to fix race condition with condition text

### QC / Coverage
- SAP #346 — Sequencing graphs
- SAP #327 — VG3 compound het node is slow (fix backport)
- SAP #395 — QC gene coverage: stop cascade deleting

### Phenotype Matching
- VG #1374 — NovaSeq X image
- Patient phenotype matching improvements: skip blank patients, remove commonly matched words

### Classification / Sync
- VG #1415 — Fix crash bug retrieving ClassificationModificationSync record
- SAP #356 — VG3 ClinVar new columns patch back
- Make djgeojson and leaflet optional (#1410)

---

## VG3 Key Characteristics (vs VG4 / master)
- Running older Django version (pre-Django 5/6 upgrade)
- Does NOT have somatic classification workflow (Horak, AMP levels)
- Does NOT have gnomAD v4 annotation
- Does NOT have T2T/CHM13v2 genome build support
- Does NOT have Structural Variant v2 improvements
- Does NOT have Candidate Search / Reanalysis tool
- Does NOT have Conservation Node
- Receives only targeted bug fixes and SAP-specific cherry-picks
- Still using VG3-era analysis pipeline and infrastructure
