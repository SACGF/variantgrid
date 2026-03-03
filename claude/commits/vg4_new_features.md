# VG4 New Features for SA Path — Medical Scientist Summary

This document summarises what is new in VG4 compared to the VG3 SA Path instance.
It focuses on features that directly affect how medical scientists classify, review, and analyse variants.

## Structural Variants (SV/CNV) Support — Major Improvements

VG3 used explicit ref/alt bases. VG4 has significantly improved handling of large variants via
handling `<DEL>`, `<DUP>` and `<INV>` with arbitrary lengths.

**What this means for medical scientists:**
- SVs are now correctly **annotated** — gene impact, population frequency (gnomAD SV), ClinVar SV overlap
- Large deletions, duplications, and inversions import and display properly
- Quick links to external databases (Varsome, gnomAD) work for SVs
- Variants from Manta, Dragen, and other SV callers import more reliably
- gnomAD SV population frequency shown on variant detail page

---

## Analysis

- **MOI (Mode of Inheritance) Node** — a dedicated analysis node to filter by inheritance pattern (autosomal dominant, recessive, X-linked, de novo), uses zygosity and gene/disease associations (GenCC)
- **Conservation Node** — filter variants by evolutionary conservation scores
- **Performance improvements** — node counts, node loading, adding tags, and creating analyses from templates are all faster
- **Prevent repeated node loading** — users can no longer accidentally load a node multiple times by repeated clicking (was causing performance issues)
- **Analysis grid tags/classifications** take into account sharing and withdrawn classifications
- **Tags shown on analysis listing page**
- **EffectNode** updated to use the new pathogenicity columns; automatically switches based on the annotation version used
- **Tissue Node** works again (updated Human Protein Atlas annotation)
- **Trio Node** Proband HET option removed and replaced with SampleNode (which is faster)
- **Population Node** — can exclude variants flagged as filtered by gnomAD; gnomAD filtering allele frequency (FAF) also available as a separate filter option; shows which gnomAD version is in use
- **PhenotypeNode** — supports MONDO disease ontology, with a warning if the ontology term doesn't map to any genes
- **Analysis grid shows patient name** alongside sample name
- **Cohort node per-sample filtering** — select individual samples to include/exclude with correct counts
- **Auto-launch analysis templates** — configure templates to auto-launch when a new sample arrives, and preview which samples would match. (In VG3 this was hardcoded for Molecular Oncology only)
- **VCF upload notification** — you are now told when uploaded VCFs are being annotated, rather than showing an empty grid
- **Analysis audit log** — records nodes being added/configured so you can track what was changed and when

---

## Annotation

VG4 includes substantially newer annotation databases:

| Database | VG3 | VG4 |
|----------|-----|-----|
| gnomAD (GRCh38) | v3 | **v4** (many new fields, hemi counts, PopMax AC) |
| gnomAD (GRCh37) | v2.1 | v2.1 + extra fields |
| dbNSFP | 4.2 | **4.5** |
| AlphaMissense | Not available | **Included** (AI-predicted pathogenicity) |
| MAVE scores | Not available | **Included** (saturation genome editing data) |
| VEP | 105 | **110** (also updated Mastermind and COSMIC) |

**New annotation columns:**
- **New dbNSFP gene columns** — GDI, GDI-Phred, P(HI), GHIS, P(rec), HIPred score, Gene indispensability score, Expression (egenetics), Expression (GNF/Atlas), BioCarta Pathway, ConsensusPathDB Pathway, KEGG IDs, KEGG Pathway, GWAS traits, GO process, GO cellular, GO molecular function, BioGRID interactions, ConsensusPathDB interactions, gnomAD prob LOF intolerant/HOM/tolerant, Essential Gene (CRISPR/CRISPR2/Gene Trap)
- **ALoFT** and **NMD Escaping Variant** — new variant-level annotation columns
- **AlphaMissense** — new missense pathogenicity predictor (Google DeepMind) — useful for BP4/PP3
- **MAVEdb** — saturation genome editing scores available for supported genes

**Gene annotation column changes:**
- OMIM/HPO terms now only use direct links (was too broad due to indirect links)
- New **MONDO terms** column
- New **MOI columns** showing gene/disease inheritance classification by confidence level (e.g. for GATA2: *Definitive: ClinGen, Strong: Genomics England PanelApp*), separately for Moderate-or-above and Supportive-or-below

**Pathogenicity Predictions — updated:**
- **New rankscore columns** — AlphaMissense, BayesDel (NoAF), CADD (raw), ClinPred, REVEL, MetaLR, VEST4
  - These are scaled 0→1 based on all non-synonymous SNVs. There is 1 score per variant (not per transcript).
  - Previously we used levels (e.g. DAMAGING, BENIGN), but those used genome-wide thresholds which vary considerably per gene.
- New pathogenicity prediction scores now **autopopulate classifications**

---

## Classification

Several improvements to the day-to-day classification experience:

- **Grouping** of similar classifications (to handle some labs that have dozens of classifications for a variant)
- **Classification Dashboard** — a new summary page showing condition text match status, lab activity, and outstanding issues (replaces the old flat "issues" list)
- **Condition text matching** uses MONDO, OMIM, and HPO to auto-suggest ontology terms when you type a condition
- **GenCC gene-disease validity** integrated — helps assess pathogenicity evidence for gene-disease relationships
- **Discordance triage improvements** — clearer interface, ability to search discordance reports, expert panel information shown
- **Pending concordance preview** — see what a concordance resolution would look like before submitting
- **`Refer` and `Seen` date fields** added to CSV export
- **Multi-lab support** — if you belong to more than one lab, a dropdown lets you choose which lab to classify under
- **Invisible character detection** — classifications with invisible Unicode characters are flagged and cleaned
- **Withdrawn classification handling** improved throughout


## Somatic Classification (Horak / AMP Framework)

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

## Candidate Search / ClinVar Reanalysis Tool

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
## Data

- **Download annotated CSV/VCF** — generate annotated output files from an analysis node for offline download on the sample/VCF page
- **Samples page** — can now filter for variant type; autocompletes fixed to search all genome builds
- **Cohorts / Pedigrees** — user is now shown alongside cohort/pedigree to make it easier to find your own

---

## Genes

- **BioCommons HGVS library** — switched from the old library, includes support for `inv` HGVS notation
- **Gene Symbol page** — rearranged and updated with new annotation fields (dbNSFP gene columns, MOI data, etc.)
- **Gene/Disease** — now uses **GenCC** (which includes ClinGen) instead of just ClinGen on gene pages and Gene Grid columns

---

## Ontology

- **Multiple ontology versions** — ontology can now have multiple versions, allowing upgrades over time without breaking existing data
- **Patients phenotype MONDO matching** — patient phenotype text is also matched against MONDO terms
- **Updated to latest versions** of HPO, OMIM, MONDO, and other ontologies
- **Ontology term page** — now shows classifications and patients linked to the ontology term

---

## Search

Variant search is substantially more informative and helpful:

- **Gene-symbol HGVS search** — search using `BRCA1:c.185del` format directly
- **HGNC ID search** — search by HGNC integer identifier
- **dbSNP lookup via API** — searches rs numbers via live API lookup (was very slow searching the large local database)
- **Much better error messages** — instead of generic failures, you now get specific messages explaining why a search didn't resolve (normalisation differences, unsupported coordinates, etc.)
- **HGVS normalisation warnings** — if your searched HGVS was normalised to a different form, a warning is shown
- **Cohort/pedigree search** — find cohorts and pedigrees by name in the search bar
- **Ontology term search** — can now search for ontology terms by name
- **Manual variant entry** — improved user experience for creating novel variants via search

---

## T2T (Telomere-to-Telomere) Genome Build

VG4 supports the new **T2T-CHM13v2** reference genome as a third genome build alongside GRCh37 and GRCh38.

**What this means for medical scientists:**
- Variants can be imported in T2T coordinates
- VEP annotation and gnomAD data available for T2T
- Liftover between T2T and GRCh38/37 supported
- Most relevant for research workflows or WGS labs adopting T2T

---

## Tagging

- **Tag Colours** — tags used to be one colour per user; they can now have multiple colours that can be shared/set globally (so all lab members can share the same colour scheme)
- **Tags for Variants grid** — the grid showing all tags per variant has been restored (was disabled for performance reasons)

---

## Variant Details

- **MANE transcripts** — transcripts are now tagged as **MANE Select** (single preferred transcript per gene) or **MANE Plus Clinical**; the MANE Select transcript is highlighted and chosen as the default on variant detail pages and in classification
- **Show all locus variants** — other variants at the same locus are now shown regardless of genotype call; also links to alleles that were once on the same row in VCF files
- **Quick gene info summary** — shows a summary of gene/disease links on the variant page (click to expand for full details)
- **Tag counts** — shows unique tag counts rather than listing all tags (which could get extremely long for popular variants). Able to click + expand if you want
- **Annotation on intergenic variants** — users can now see limited annotation on intergenic variants
- **Tooltips** — added tooltips to show help about columns on gene/variant pages
- **Nearby variants** — now shows tags and classification summary counts, and works across genome builds
- **Structural Variant page** — page appears differently for SVs (unavailable annotations hidden, notes field added)
- **ClinVar** - load more detailed, up to date info from ClinVar website

---

## Variant Grids

- **Column labels** fixed and standardised; new links added for ontology terms, PubMed, and dbSNP
- **Consistent styling** — variant pages (gene symbol, all variants, tags, etc.) now styled and linked the same way as the analysis grid

---

## VCF Upload

- **Undefined filter warning** — if VCF filters are not defined in the VCF header, a warning is now given rather than failing to import the VCF
- **MT (chrM) variants** — various fixes so mitochondrial variants are no longer skipped or mishandled during import
- **Custom VCF INFO/FORMAT fields** — custom annotation fields from VCF files can now be stored and used in analysis (read/stored; display in grids coming later)
