# Paywalled / hard-to-retrieve articles

These are references relevant to the VariantGrid paper that I could **not** download as
open-access full text. For each I have read the abstract / publisher landing page / GitHub
or documentation. Please retrieve the full-text PDF (institutional access) and drop it into
`paper/literature/pdfs/` using the suggested filename, then re-run `python3 convert.py`.

> Status legend: **NEEDED** = please retrieve; **NICE-TO-HAVE** = abstract is probably
> enough for our purposes but full text would help the comparison/discussion.

---

## Standards & guidelines (core citations — NEEDED)

### ACMG/AMP variant interpretation guidelines
- Richards S, et al. "Standards and guidelines for the interpretation of sequence variants:
  a joint consensus recommendation of the American College of Medical Genetics and Genomics
  and the Association for Molecular Pathology." *Genetics in Medicine* 17.5 (2015): 405–424.
- DOI: 10.1038/gim.2015.30
- Why: the foundational classification framework VariantGrid implements. Must cite.
- Suggested filename: `acmg_richards_2015.pdf`

### ClinGen — Clinical Genome Resource
- Rehm HL, et al. "ClinGen — the Clinical Genome Resource." *New England Journal of Medicine*
  372.23 (2015): 2235–2242.
- DOI: 10.1056/NEJMsr1406261  (NEJM — paywalled)
- Why: the curation/standards body behind discordance resolution and variant curation expert
  panels; context for Shariant/discordance work.
- Suggested filename: `clingen_rehm_2015.pdf`

### Sherloc / refinements of ACMG (NICE-TO-HAVE)
- Nykamp K, et al. "Sherloc: a comprehensive refinement of the ACMG–AMP variant
  classification criteria." *Genetics in Medicine* 19.10 (2017): 1105–1117.
- DOI: 10.1038/gim.2017.37
- Why: relevant to the flexible evidence-key / criteria-strength system.
- Suggested filename: `sherloc_nykamp_2017.pdf`

---

## Comparator platforms with no open-access primary paper

### seqr — final journal version (preprint IS open access; see pdfs/seqr_2022)
- Pais LS, et al. "seqr: A web-based analysis and collaboration tool for rare disease
  genomics." *Human Mutation* 43.6 (2022): 698–707.
- DOI: 10.1002/humu.24366  (Wiley — paywalled; we have the medRxiv preprint already)
- Why: nearest open-source comparator. Preprint is sufficient for comparison; final version
  NICE-TO-HAVE only.
- Suggested filename: `seqr_hummut_2022.pdf`

### Commercial platforms — no peer-reviewed "tool paper" exists
These are cited from vendor documentation / white papers, not journal articles. No retrieval
needed, listed so we remember to cite them correctly as software/URLs:
- **Agilent Alissa Interpret** (formerly Cartagenia Bench Lab) — vendor product page only.
- **Golden Helix VarSeq** — vendor product page / application notes only.
- **Genoox Franklin** — vendor product page; appears in the 7-platform comparison we have.
- **Illumina Emedgene** — vendor product page; appears in the 7-platform comparison.
- **QIAGEN Clinical Insight (QCI) Interpret** — vendor product page; in the comparison.
- **VarSome Clinical (Saphetor)** — has a primary paper (Kopanos 2019, Bioinformatics);
  see open-access list, retrieved separately if possible.

---

## Reclassification / reanalysis clinical-impact studies (NICE-TO-HAVE)

These support the "genetic results have a long lifetime; reanalysis pays off" argument in the
introduction/discussion. Abstracts are captured; full text strengthens the numbers we quote.

### Diagnostic yield of exome reanalysis over time
- *Genetics in Medicine* (2026), S1098-3600(26)00867-1.
- Why: timing/yield data for periodic reanalysis.
- Suggested filename: `reanalysis_yield_over_time_gim.pdf`

### Cardiovascular / arrhythmia variant reclassification impact
- VanDyke, et al. *Journal of Genetic Counseling* 30 (2021). DOI: 10.1002/jgc4.1336
- Inherited arrhythmia clinic reclassification frequency — *Heart Rhythm* (2024).
- Why: clinical-management-change statistics for the reclassification motivation.
- Suggested filenames: `reclass_cardio_vandyke_2021.pdf`, `reclass_arrhythmia_2024.pdf`

---

## Have as web-extract only — full PDF still wanted (NICE-TO-HAVE)

These could NOT be auto-downloaded as PDFs (publisher Cloudflare / EuropePMC render failures),
but I captured a structured text extract of each in `paper/literature/text/*_EXTRACT.txt`. The
extracts are good enough for the comparison; retrieve the full PDF only if we end up quoting
methods/figures in detail. **VarFish is the priority** of this group (closest comparator).

| Paper | Extract file | Full-PDF source to fetch |
|-------|--------------|--------------------------|
| **VarFish** (Holtgrewe 2020, NAR W1) | `varfish_2020_EXTRACT.txt` | doi:10.1093/nar/gkaa241 / PMC7319464 |
| seqr (Pais 2022, Hum Mutat) | `seqr_2022_EXTRACT.txt` | doi:10.1002/humu.24366 (preprint medRxiv 2021.10.27.21265326 is OA) |
| iVar (2021, *Genes*) | `ivar_2021_EXTRACT.txt` | PMC8001268 (MDPI, OA — should be retrievable) |
| Matchmaker Exchange (Philippakis 2015) | `matchmaker_decipher_2022_EXTRACT.txt` | PMC4610002 |
| DECIPHER (Foreman 2022, Hum Mutat) | `matchmaker_decipher_2022_EXTRACT.txt` | PMC9303633 |
| 7-platform comparison (2025) | `platform_compare_2025_EXTRACT.txt` | PMC11949535 (PDF downloaded but did not text-convert — `pdfs/platform_compare_2025.pdf` is on disk) |

## Successfully retrieved as full text (for reference — no action needed)

In `paper/literature/text/` as full PDF→text conversions: shariant_2022, vep_2016, gemini_2013,
vcfminer_2016, opencravat_2019, clinvar_2018, exomiser_2015, exomiser_reinterp_2024,
exome_reanalysis_review_2021, cancer_reclass_2020.
</content>
