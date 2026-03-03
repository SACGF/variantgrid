# VG4 New Features for SA Path — Medical Scientist Summary

This document summarises what is new in VG4 compared to the VG3 SA Path instance.
It focuses on features that directly affect how medical scientists classify, review, and analyse variants.

---

## 1. Somatic Classification (Horak / AMP Framework)

VG4 has full support for classifying somatic/oncology variants using the **Horak (AMP/VICC) framework** alongside the existing germline ACMG system.

**What this means for medical scientists:**
- Classifications now have an **Allele Origin** (Germline / Somatic / Unknown) — this is mandatory
- Somatic variants use **Oncogenic / Likely Oncogenic / VUS / Likely Benign / Benign** classification
- AMP Levels (Tier I–IV) are displayed throughout the system — in classification lists, search results, and the variant page
- The classification form adapts based on allele origin: somatic fields appear for somatic records, germline fields for germline
- Assertion method can now be both ACMG and Horak simultaneously (for labs classifying both)
- ClinVar oncogenic data is imported and shown alongside germline data
- Classification lists clearly separate germline and somatic columns

---

## 2. Candidate Search / ClinVar Reanalysis Tool

A brand-new SA Path workflow to **find variants in your samples that may need reclassification** based on updated ClinVar evidence or new publications.

**What this means for medical scientists:**
- Search across your sample variants filtered by gene/phenotype/zygosity
- The tool highlights variants where ClinVar has changed since your last classification
- Shows your existing classifications alongside current ClinVar data for comparison
- You can **classify directly from the search results** without needing to go back to the full analysis
- Filter by individual sample, VCF project, or phenotype/ontology term
- Shows how long ago the variant was last classified
- Medical scientists can change workflow state directly from the reanalysis pages (not just admins)

---

## 3. Structural Variant (SV/CNV) Support — Major Improvements

VG3 used explicit ref/alt bases. VG4 has significantly improved handling of large variants via
handling <DEL>, <DUP> and <INV> with arbitrary lengths

**What this means for medical scientists:**
- SVs are now correctly **annotated** — gene impact, population frequency (gnomAD SV), ClinVar SV overlap
- Large deletions, duplications, and inversions import and display properly
- Quick links to external databases (Varsome, gnomAD) work for SVs
- Variants from Manta, Dragen, and other SV callers import more reliably
- ClinVar matching for SVs uses **overlap percentage** rather than exact coordinates (more clinically meaningful)
- gnomAD SV population frequency shown on variant detail page

---

## 4. Updated Annotation Data

VG4 includes substantially newer annotation databases:

| Database | VG3           | VG4 |
|----------|---------------|-----|
| gnomAD (GRCh38) | v3            | **v4** (many new fields, hemi counts, PopMax AC) |
| gnomAD (GRCh37) | v2.1          | v2.1 + extra fields |
| dbNSFP | 4.2           | **4.5** |
| AlphaMissense | Not available | **Included** (AI-predicted pathogenicity) |
| MAVE scores | Not available | **Included** (saturation genome editing data) |
| VEP | 105           | **110** |
| ClinVar | VCF only      | **XML** (richer data, more variants) |

**What this means for medical scientists:**
- More accurate population frequency data for rare variant filtering (gnomAD v4)
- **AlphaMissense** provides a new missense pathogenicity predictor (Google DeepMind, 2023) — useful for BP4/PP3
- MAVE saturation genome editing scores available for supported genes
- ClinVar data is more complete (includes multi-allelic records and more structured condition data)

---

## 5. MANE Transcript Support

VG4 uses the **MANE (Matched Annotation from NCBI and EMBL-EBI)** transcript standard.

**What this means for medical scientists:**
- Transcripts are now tagged as **MANE Select** (single preferred transcript per gene) or **MANE Plus Clinical**
- The classification form highlights the MANE Select transcript to guide transcript choice
- Hotspot graphs can be viewed on the MANE Select transcript
- Consistent transcript selection across the system reduces discordance between labs

---

## 6. T2T (Telomere-to-Telomere) Genome Build

VG4 supports the new **T2T-CHM13v2** reference genome as a third genome build alongside GRCh37 and GRCh38.

**What this means for medical scientists:**
- Variants can be imported in T2T coordinates
- VEP annotation and gnomAD data available for T2T
- Liftover between T2T and GRCh38/37 supported
- Most relevant for research workflows or WGS labs adopting T2T

---

## 7. Improved Classification Workflow

Several improvements to the day-to-day classification experience:

- **Classification Dashboard** — a new summary page showing condition text match status, lab activity, and outstanding issues (replaces the old flat "issues" list)
- **Condition text matching** uses MONDO, OMIM, and HPO to auto-suggest ontology terms when you type a condition
- **GenCC gene-disease validity** integrated — helps assess pathogenicity evidence for gene-disease relationships
- **Discordance triage improvements** — clearer interface, ability to search discordance reports, expert panel information shown
- **Pending concordance preview** — see what a concordance resolution would look like before submitting
- **`Refer` and `Seen` date fields** added to CSV export
- **Multi-lab support** — if you belong to more than one lab, a dropdown lets you choose which lab to classify under
- **Invisible character detection** — classifications with invisible Unicode characters are flagged and cleaned
- **Withdrawn classification handling** improved throughout

---

## 8. Analysis Pipeline Improvements

Improvements to the variant filtering analysis:

- **MOI (Mode of Inheritance) Node** — a dedicated analysis node to filter by inheritance pattern (autosomal dominant, recessive, X-linked, de novo)
- **Conservation Node** — filter variants by evolutionary conservation scores
- **MONDO phenotype matching** — PhenoNode now matches against MONDO disease ontology in addition to OMIM/HPO
- **Analysis grid now shows patient name** alongside sample name (previously only sample ID)
- **Cohort node per-sample filtering** — select individual samples to include/exclude with correct counts
- **gnomAD filtering allele frequency (FAF)** available as a filter option in Population Node
- **Auto-launch analysis templates** — you can now configure templates to auto-launch when a new sample arrives, and preview which samples would match the configuration. This was added to VG3 but is hardcoded just for Molecular Oncology.
- **VCF upload notification** — you are now told when uploaded VCFs are being annotated, rather than showing an empty grid

---

## 9. Search Improvements

Variant search is substantially more informative and helpful:

- **Gene-symbol HGVS search** — search using `BRCA1:c.185del` format directly
- **HGNC ID search** — search by HGNC integer identifier
- **dbSNP lookup via API** — searches rs numbers via live API lookup as it was very slow searching large local DB
- **Much better error messages** — instead of generic failures, you now get specific messages explaining why a search didn't resolve (normalization differences, unsupported genome build coordinates, etc.)
- **HGVS normalisation warnings** — if your searched HGVS was normalised to a different form, a warning is shown
- **Cohort/pedigree search** — you can now find cohorts and pedigrees by name in the search bar

---

## 10. Other

- **MT chrM mito variants** — Various fixes related to MT variants
- **Custom VCF INFO/FORMAT fields** — VCF custom annotation fields can now be stored and used in analysis. These are read/stored but not displayed yet.
- **Annotated VCF/CSV download** — generate annotated output files from an analysis node for offline download on sample/VCF page
