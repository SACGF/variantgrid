# Master Branch — Feature Development Summary (VG4)

**Branch:** `master`
**Diverged from VG3 at:** 2021-06-22 (commit `62035452d`)
**Unique commits (not in vg3_sapath_prod):** 5,599
**Date range:** 2021-06-22 → 2026-02-27

This represents the VG4 development track. The scope below covers features developed since the VG3/master divergence.

---

## 1. Somatic Classification (Horak / AMP Framework)

**Period:** 2023-2024 (major work Q1-Q3 2024)

- Full AMP/VICC (Horak) somatic classification framework alongside germline ACMG
- `allele_origin` field (germline / somatic / unknown) as first-class concept
- AMP Level displayed throughout classification grids, overlaps, and quick views
- Oncogenic ClinVar data imported (`ONC` entries from ClinVar XML)
- Somatic Clinical Significance separate from germline Clin Sig in all grids
- Classification/ClinSig displayed as separate pills in overlaps page
- Somatic evidence keys (testing_context, tumour type, etc.)
- Assertion method now multi-select (supports ACMG + Horak simultaneously)
- ClinVar export updated to support somatic submissions
- "No-Data" shown appropriately based on allele origin
- Hide germline-specific columns in somatic view and vice versa
- SAP #369 — Candidate Search / Reanalysis Tool (full workflow, see below)

---

## 2. Structural Variants (SV) — Major Improvements

**Period:** 2022-2025

- SV annotation v2: gnomAD SV overlap, ClinVar SV overlap
- BCFTools-based normalisation replacing VT (handles symbolic alts properly)
- Symbolic alt support: `<DEL>`, `<DUP>`, `<INV>`, `<DUP:TANDEM>`
- SVLEN handling improvements (Manta, Dragen, bcftools edge cases)
- Variants >1kb stored as symbolic alts for efficiency
- VEP annotation of SVs (with size limits for very large variants)
- SV read length stats / graph
- Quick links working for symbolic alts (Varsome, gnomAD)
- Liftover of symbolic variants via BCFTools liftover
- ClinVar SV overlap matching by % rather than exact
- gnomAD SV population frequency on variant details
- Search support for symbolic variant HGVS

---

## 3. T2T / CHM13v2 Genome Build Support

**Period:** late 2024

- T2T-CHM13v2 as a supported genome build
- Annotation pipeline support (VEP, gnomAD liftover, RepeatMasker)
- VCF import with T2T common variant filters
- Allele liftover between T2T and GRCh37/38
- ClinVar import for T2T
- Classification grid support for T2T records

---

## 4. Annotation — gnomAD v4, AlphaMissense, dbNSFP 4.5, VEP 110

**Period:** 2023-2024

- gnomAD v4 annotation (GRCh38): new fields, hemi counts, PopMax AC
- gnomAD v2 extra fields for GRCh37
- AlphaMissense pathogenicity scores integrated via dbNSFP and VEP plugin
- dbNSFP 4.5 upgrade (raw pathogenicity scores alongside rank scores)
- VEP 110 new fields added to columns
- MAVE (MaveDB saturation genome editing) scores annotation (#850)
- Conservation Node for analysis (new analysis filter node)
- COSMIC version now tracked in annotation version
- Annotation run management improvements (kick off from page, subdivide stuck runs)
- ClinVar XML-based import (richer data, more variants vs VCF-only)

---

## 5. Candidate Search / ClinVar Reanalysis Tool (SAP #369)

**Period:** Nov-Dec 2025

A major new workflow for SA Path:
- Search for variants in your samples that have updated evidence/classifications
- Filter by sample, VCF/project, zygosity, ontology/phenotype
- ClinVar reanalysis: find your samples' variants with new ClinVar data
- View existing classifications on those variants
- Action button to create new classification directly from search results
- Users can change workflow state from the reanalysis pages
- Shows date since last classified, links to analysis

---

## 6. MANE Transcript Support

**Period:** 2022-2023

- MANE Select and MANE Plus Clinical transcript tags
- cdot-based transcript data (replacing direct Ensembl/RefSeq loading)
- MANE transcript selection in classification form
- Hotspot graph transcript picker supports MANE
- MANE tags shown on transcript pages and classification form
- Import cdot latest management command

---

## 7. Classification Improvements

**Period:** 2021-2025

- Dashboard page (replaced old "issues" page) with condition text match counts
- Condition text matching to ontology (MONDO, OMIM, HPO) auto-matching
- GenCC gene-disease validity integrated into condition text matching
- Discordance report triage workflow
- VUS (a/b/c) VUS overlap classification groupings
- Expert panels shown on discordance page
- Discordance report search
- Pending concordance preview
- `allele_origin` mandatory for all records going forward
- Classification select-lab dropdown when user belongs to multiple labs
- Bulk change sharing levels
- Classification grouping popup showing all classifications for a grouping
- Diff across allele origins on allele page
- AMP Level shown throughout
- Withdrawn classification handling improvements
- Remove non-ASCII/invisible characters from classification text
- `refer` and `seen` date fields in CSV export
- ClinVar SCV search (#3556)
- Lab compare export improvements (#3539)

---

## 8. Search Improvements

**Period:** 2022-2025

- HGVS search completely rewritten: better error messages, less cryptic failures
- Gene-symbol HGVS search (e.g. `BRCA1:c.1234A>T`)
- HGNC integer search + link to gene symbol page
- dbSNP search via API lookup
- Search for cohort/pedigree
- MT HGVS search across genome builds
- hg18/patch version coordinates give informative error
- HGVS normalization warning shown in search results
- Search results grouped/collated to avoid duplicates
- Symbolic variant search improvements
- Incomplete HGVS locus handling (e.g. `NM_015120.4:c.4166`)
- `c.` without transcript gives informative error

---

## 9. Analysis Pipeline Improvements

**Period:** 2021-2025

- MOI (Mode of Inheritance) analysis node (#548) — new node type
- Conservation Node — new node type (filters on conservation scores)
- MONDO added to Phenotype node
- Analysis grid CSV download from output node (#1398)
- Cohort node: per-sample zygosity warnings, correct counts on individual-sample filter
- Analysis grid shows patient name alongside sample name (#1318)
- Population node: gnomAD help, filtering allele frequency option
- Analysis templates: view how samples would match auto-run templates (#1428)
- Analysis auto-launch: configuration and auto-launch management (#1384)
- Candidate Search / Reanalysis (see §5)
- VCF upload: user notified that annotation is in progress (#1269)
- Custom cohort: handles multiple cohorts, copy common genotypes (#3704)
- Annotated VCF / CSV generation for download (#1171)

---

## 10. Sequencing Automation (SeqAuto) — VG4 Additions

**Period:** 2022-2025

- SeqAuto API v2 (dataclasses-based, gene coverage, QC)
- Single-sample VCF support via SeqAuto API
- NovaSeq X support (sequencer creation via API)
- Gene coverage reload stability improvements
- Sample GOI duplication fix (API + scan)
- Custom VCF INFO/FORMAT fields support (#1240)
- Reload VCF in-place capability (#76)
- Optimised VCF sample stats (#1250)

---

## 11. Django / Infrastructure Upgrades

**Period:** 2022-2026

- Django 4 (Jan 2022)
- Django 5 (Dec 2025)
- Django 6 (current, 2026)
- Python 3.11 → 3.12 support
- JQuery Layout replaced with modern CSS (#978)
- JQGrid-to-DataTables conversion ongoing (custom columns, ClinVar export done)
- Celery settings modernised
- Dependencies upgraded throughout (pyVCF2, PyPDF, django-avatar, captcha)

---

## 12. Patients / Phenotype

**Period:** 2021-2025

- MONDO added to phenotype matching
- MDS (Multi-Dimensional Scaling) phenotype matching improvement (#1203)
- Patient phenotype matching: skip blanks, speed improvements
- Patient page shows samples on cohort/trio pages (#538)
- Patient records import improvements (#1086)
- Patient import handles new fields, deleted samples (#1327)
- Prenatal age field (#1449 partial)

---

## 13. Liftover / Genome Build

**Period:** 2022-2025

- BCFTools liftover (replacing Picard) for classification imports
- Symbolic variant liftover handling
- AlleleLiftover failure events tracked
- T2T liftover support (see §3)
- Proper handling of genome build patch versions (e.g. GRCh38.p14)
- hg19 vs GRCh37 handling improvements

---

## 14. ClinVar Export

**Period:** 2022-2025

- ClinVar XML import (richer data)
- Orphanet → ORPHA format fix for ClinVar submission
- ClinVar export supports both germline and somatic (#3593)
- Orphaned ClinVar export record cleanup (#3165)
- Better display of withdrawn ClinVar exports (#3341)
- ClinVar batch health checks (#3554)
- Expert panel data on discordance page (from ClinVar XML)
- NC_ transcript support for ClinVar submission

---

## 15. Minor Quality-of-Life

- Emojis in EventLog (#1098)
- Tag colours more user-friendly (#2439)
- Tags shown in nearby variants summary (#137)
- Tags shown on analysis list page
- Nicer Gene ID page formatting (#803)
- Site messages support Markdown (#1378)
- VCF filter tag display improvements
- Read-only icon on listing pages (#968 partial)
- Annotation version shows cdot version (#958)
- Server status page improvements (worker status, Celery)
- Pedigree file deletion
- Gene wiki sortable by gene symbol
