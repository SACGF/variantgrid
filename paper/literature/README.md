# VariantGrid paper — literature folder

Literature search supporting the VariantGrid platform paper. Generated 2026-06-12.

## Contents
- `pdfs/` — downloaded open-access PDFs.
- `text/` — readable text. `*.txt` = full PDF→text conversions; `*_EXTRACT.txt` = structured
  web extracts for papers whose PDF could not be auto-downloaded (see `paywalled_articles.md`).
- `convert.py` — re-run `python3 convert.py` after dropping new PDFs into `pdfs/`.
- `comparison.md` — **the main analysis**: landscape, comparator-by-comparator table, where
  VariantGrid is novel, where claims need softening.
- `feature_comparison.csv` — **feature matrix** (~65 features × 9 tools). VG column source-verified
  with file references; comparator columns from the literature (`?` = unverified).
- `feature_comparison_notes.md` — how to read the CSV, the STRONG/PARITY/BEHIND groupings, and a
  proposed ~15-row headline table for the paper.
- `paywalled_articles.md` — references to retrieve manually (institutional access), incl. the
  must-cite standards (ACMG Richards 2015, ClinGen Rehm 2015).
- `../journal_targets.md` — ranked candidate journals + recommendation.
- `../2018_paper.txt` — the old (out-of-date) draft.

## Sources retrieved (15 covered)

**Closest open-source comparators:** VarFish (extract), seqr (extract), GEMINI, VCF-Miner,
OpenCRAVAT. **Annotation:** VEP, (ANNOVAR — paywalled list). **Sharing/standards ecosystem:**
Shariant (VariantGrid's own national deployment), ClinVar, Matchmaker Exchange + DECIPHER
(extract). **Re-analysis/re-classification:** Talos (Centre for Population Genomics 2025 — key automated-
reanalysis comparator, complementary not rival), iVar (extract — direct prior art),
exome-reanalysis review, Exomiser + Exomiser-reinterpretation, cancer variant reclassification.
**Commercial landscape:** 7-platform comparison (extract).

## One-line takeaway

VarFish (PostgreSQL + Django, on-prem, ACMG, re-analysis) is the single closest comparator and
the must-have head-to-head row. VariantGrid's defensible novelty is the **combination** of
(1) live node-graph filtering, (2) annotation-version partitioning enabling constant-time
multi-version analysis + automated re-analysis by version-diffing, (3) SQL-native multi-sample
array packing, and (4) the end-to-end multi-lab classification-sharing/discordance/ClinVar
ecosystem proven at national scale by Shariant. Recommended venue: **NAR Web Server issue**
(where VarFish published), backup **Bioinformatics**.
