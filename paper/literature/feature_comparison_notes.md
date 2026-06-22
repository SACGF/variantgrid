# Feature comparison — notes & how to read it

Companion to `feature_comparison.csv`.

## Method & confidence
- **VariantGrid column is source-verified** — every Y has a file/class reference in the
  `vg_evidence_or_note` column, from a code survey of `snpdb/ upload/ annotation/ analysis/
  classification/ sync/ flags/ ontology/ variantopedia/` (2026-06).
- **Comparator columns are from the literature** in `paper/literature/` (VarFish, seqr, Scout,
  GEMINI, VCF-Miner, OpenCGA, Alissa, Franklin/VarSome). Cells marked `?` are unverified — I did
  not have authoritative text. **Verify the `?` and `P` comparator cells before publication**
  (vendor docs / tool papers), since a feature table is the part reviewers scrutinise hardest.
- `vg_differentiator` summarises how VG compares per feature:
  `STRONG` = unique/near-unique to VG · `MOD` = VG does it more/better than most ·
  `PARITY` = common, VG has it like others · `BEHIND` = VG weaker than some comparators.

## The STRONG cluster (VG's defensible novelty), grouped

**1. Database / scale backbone**
- Many VCFs → one cross-build `Allele` record; multi-sample genotype packing for SQL-level
  filtering (extends GEMINI into Postgres arrays); sample-level "who carries this" traceback;
  rare/common AF-based partition that skips the common partition for rare queries; data-archive
  lifecycle for long-term retention without active-query cost.

**2. Annotation versioning**
- Multiple annotation versions stored at once; **partition-per-version** for constant-time
  analysis; **version diffing with transition matrices** driving re-analysis; per-analysis/
  per-classification version pinning. No comparator does this set.

**3. Live node-graph analysis**
- Interactive DAG filtering with live editing/immediate feedback; ~20 node types; Venn/set ops;
  auto-launch on import; node count badges. Alissa is node-based but **batch**; everyone else is
  form/query based.
- **Alissa Interpret (the node-based pioneer, ex-Cartagenia Bench) has been END-OF-LIFED** —
  Agilent announced discontinuation ~2023, retired end of 2024; the vendor migration path (QCI
  Interpret) is not node-based. So VariantGrid is the **only actively-maintained node-graph
  variant analysis platform**, and the only open-source/live one. (Cite as Agilent/QIAGEN vendor
  EOL notice + accessed date — not a formal publication.) This is a strong, clean claim; lead
  with "only actively-maintained / only open-source live," not a bare "first."

**4. Multi-build linking + classification-sharing ecosystem**
- Cross-build linking of variants AND classifications/tags; lab-configurable evidence keys;
  multi-lab sharing with share levels; automated inter-lab discordance detection + triage +
  per-condition clinical grouping; ClinVar pipeline; two-way Shariant sync; condition→ontology
  matching; internal variant matching.

## Where VG is at PARITY (don't claim novelty)
Open source, REST API, ACMG workflow, trio/cohort/pedigree, gene/PanelApp panels, in-silico
predictors, HGVS search, variant tagging, object permissions. These are table stakes — include
them in the table for completeness, but the paper's claims should rest on the STRONG cluster.

## Where VG is BEHIND (be honest; pre-empt reviewers)
- **Phenotype-based ranking:** VarFish/seqr/Scout/Exomiser/commercial rank variants by
  phenotype similarity; VG has a `PhenotypeNode` filter but **not a ranking/prioritisation
  engine**. If reviewers expect Exomiser-style ranking, acknowledge it (and that VG can consume
  external predictors) rather than overclaim.
- **Matchmaker Exchange:** seqr and Scout integrate MME for cross-institution gene discovery;
  VG does **internal** (within-instance) variant matching only. Frame as complementary, not a
  gap to hide. (Possible future-work line.)
- **Docker/turnkey deploy:** comparator parity at best — VG stack is heavy; Docker image is
  planned, not shipped. Don't claim ease.

## Suggested table for the paper
The CSV has ~65 features — too many for a main-text table. For the paper, collapse to a
**~15-row headline table** drawn from the STRONG cluster + a few PARITY rows for context, columns
= VariantGrid · VarFish · seqr · Scout · OpenCGA · Alissa (the closest functional peers). Keep the
full CSV as supplementary. VarFish is the must-include column (closest peer).
