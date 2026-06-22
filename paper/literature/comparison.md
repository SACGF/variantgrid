# VariantGrid — Literature Landscape & Positioning

_Prepared as background for the VariantGrid platform paper. Summarises the comparable tools
retrieved into `paper/literature/`, what each does, and where VariantGrid is genuinely novel vs.
where it overlaps. Sources are in `pdfs/` (full PDFs) and `text/` (converted text or web
extracts where the PDF could not be auto-downloaded — see `paywalled_articles.md`)._

## 1. What VariantGrid is (for positioning)

A self-hostable Django/PostgreSQL web platform that spans the **whole** variant lifecycle:
VCF upload + normalisation → managed VEP annotation → interactive node-graph filtering/analysis
→ ACMG/AMP classification → multi-lab classification sharing, discordance detection, and ClinVar
submission. Deployed as variantgrid.com (research), SA Pathology (clinical), RUNX1db (disease
database), and **Shariant** (Australian national classification-sharing network).

Four claimed novelties (from the 2018 draft, still the core of the story):
1. **Live, interactive node-graph (DAG) variant filtering** with immediate feedback.
2. **Partition-per-annotation-version storage** → constant analysis performance regardless of
   how many annotation versions are stored, enabling **automated re-analysis/re-classification**
   by diffing annotation versions.
3. **Gemini-style multi-sample packing extended to Postgres array columns**, so filtering across
   thousands of samples stays in SQL (the DB optimiser can reorder it).
4. **An end-to-end classification + multi-lab sharing + discordance + ClinVar ecosystem**
   (Shariant), not just an analysis front-end.

A fifth, database-side differentiator (raised by D. Lawrence, and present in the 2018 draft but
under-sold):
5. **Continuous operational ingestion of routine VCFs, each normalised and linked to ONE
   cross-build `Allele` record** — building an always-growing institution-internal population
   resource with **sample/patient-level traceback** ("who else carries this allele"), i.e.
   automatic *within-institution* matchmaking + an internal allele-frequency annotation that
   catches what gnomAD cannot (pipeline artifacts, local founder variants).

Most tools below own one or two of these slices; almost none own all five, and the combination
plus the annotation-versioning design is the defensible contribution.

### Note on differentiator #5 — be precise vs comparators
- **seqr** is organised around per-project / per-family **joint-called callsets**, not a rolling
  institution-wide allele record; it annotates with *external* frequencies and does not surface
  "which of our samples carry this allele." → clean differentiator vs seqr.
- **VarFish** DOES let an institution "build a database of variants identified at the user's
  institution," so internal frequency itself is **not** unique. The finer, defensible distinction
  is: (a) always-on ingestion into one normalised cross-build allele record (not a per-analysis
  callset), (b) **sample-level traceback**, not just an aggregate count, and (c) it scales via the
  SQL-array packing (2018 draft: 668 exomes → 2.4B genotypes vs 46M variants, ~50×) — the same
  evidence ties #3 and #5 into one argument.
- Complementary to **Matchmaker Exchange** (cross-institution discovery): VariantGrid does the
  *within*-institution version automatically as a by-product of ingestion.

---

## 2. Comparator landscape (grouped)

### A. Open-source web variant analysis/interpretation platforms (closest competitors)

| Tool | Stack | Scope | Closest overlap | Key difference from VariantGrid |
|------|-------|-------|-----------------|----------------------------------|
| **VarFish** (Holtgrewe 2020, NAR) | **PostgreSQL + Django** | filter/prioritise/annotate, ACMG user classification, on-prem, case re-analysis | **Very high** — nearly identical stack & feature surface | Form-based SQL query builder (not a live node-graph); star schema (not partition-per-annotation-version); single-institution (no inter-lab discordance/sharing); re-analysis = re-run vs updated background DBs, not versioned diffs |
| **seqr** (Pais 2022, Hum Mutat) | Python, GCP/Kubernetes | family-based rare disease analysis + collaboration + Matchmaker | High (analysis + collaboration) | Cloud-native (vs self-hostable Postgres); research-only (no PHI); form-based search; no annotation-version snapshots; strong on matchmaking, weak on lab-grade classification/discordance |
| **Scout** (Clinical Genomics, Stockholm) | Python/Flask + MongoDB | clinical VCF visualiser/browser, case collaboration, Matchmaker | Medium (case review UI) | Document DB; visualiser/browser emphasis, not a scalable shared variant warehouse or node-graph; no annotation-versioning |
| **GEMINI** (Paila 2013, PLoS Comp Biol) | SQLite, CLI/Python | query genetic variation + annotations; multi-sample packing | Medium (the packing idea VariantGrid extends) | CLI library, not a web app/clinical workflow; SQLite single-file; packed data extracted in Python (VariantGrid's key improvement: keep it in SQL) |
| **VCF-Miner** (Hart 2016, Brief Bioinform) | MongoDB, GUI | stepwise GUI mining of VCF + annotations | Low–medium (filtering UI) | User-supplied annotations (no managed annotation); no classification/sharing; no versioning |
| **OpenCGA / OpenCB** (Zetta Genomics) | Hadoop/Spark/MongoDB/Solr | big-data variant storage + clinical interpretation at 100k-genome scale | Medium (scalable storage) | Heavyweight big-data infra; storage/interpretation engine rather than an interactive analyst-facing workflow; not aimed at small-lab self-hosting |

### B. Annotation engines (VariantGrid consumes/manages these, not a competitor)

| Tool | Note |
|------|------|
| **Ensembl VEP** (McLaren 2016, Genome Biol) | VariantGrid's annotation engine. Cite as dependency. |
| **ANNOVAR** (Wang 2010, NAR) | Classic annotation tool; context for "managed vs user annotation" discussion. |
| **OpenCRAVAT** (Pagel 2019/2020) | Modular meta-annotator; alternative philosophy (user-composed modules) vs VariantGrid's managed, versioned annotation. |

### C. Classification / evidence-sharing ecosystem (where Shariant sits)

| Resource | What it shares | Relationship to VariantGrid/Shariant |
|----------|----------------|--------------------------------------|
| **Shariant** (Stark/Spurdle 2022, AJHG) | **Structured classification evidence between Australian labs**, discordance resolution, ClinVar streamlining | This **is** VariantGrid deployed as a national network — the platform paper should cite Shariant as the flagship deployment/outcome, and avoid double-counting it as a "competitor." |
| **ClinVar** (Landrum 2018, NAR) | Aggregated variant interpretations | VariantGrid submits to it; positions VariantGrid as upstream structured-evidence layer. |
| **ClinGen** (Rehm 2015, NEJM) | Curation standards, VCEPs | Standards body; context for discordance/criteria work. |
| **DECIPHER** (Foreman 2022, Hum Mutat) | Anonymised phenotype-linked variant data, proportionate consent | Sharing-model parallel (share levels), but discovery/matchmaking focus, not lab classification concordance. |
| **Matchmaker Exchange** (Philippakis 2015, Hum Mutat) | Federated phenotype+genotype matching for gene discovery | Complementary sharing goal (discovery), not classification concordance. |

### D. Commercial platforms (cite as landscape; mostly no primary paper)

Agilent **Alissa Interpret** (ex-Cartagenia Bench — the first node-based variant filter UI, but
**batch** execution; VariantGrid's live editing is the differentiator. **NOW END-OF-LIFED** —
Agilent retired it (announced ~2023, end 2024); users migrated to QCI Interpret, which is not
node-based — so VariantGrid is the only actively-maintained node-graph tool), Golden Helix **VarSeq**,
Genoox **Franklin**, Illumina **Emedgene**, **QIAGEN QCI**, **VarSome Clinical** (Kopanos 2019).
The 7-platform comparison (PMC11949535, 2025) benchmarks several of these on classification &
prioritisation — useful to cite for (a) the crowded commercial field and (b) that
classification concordance is the field's key evaluation axis, which Shariant addresses head-on.
None are open-source/self-hostable → reinforces VariantGrid's data-sovereignty positioning.

### E. Re-analysis / re-classification (motivation AND the key comparator: Talos)

- **Talos** (Centre for Population Genomics, medRxiv 2025; `talos_2025.txt`) — current
  state-of-the-art **automated re-analysis**. Australian, open-source, **batch pipeline**
  (Nextflow/Hail/Docker, cohort-scale, monthly ClinVar/PanelApp cycles, JSON output); 5.2% extra
  yield on 4,735 undiagnosed; compares itself to Exomiser. **Complementary, not a competitor:**
  Talos is a standalone flagging pipeline with no UI/persistent queryable DB/interactive analysis;
  VariantGrid is a platform where reanalysis is integrated (annotation-version diffing inside a DB
  that holds the samples, classifications, tags, analyses) plus live node-graph exploration. They
  could interoperate. **Do NOT try to out-yield Talos** — position VG as *infrastructure for
  reanalysis*, not a rival yield engine. (See `../strategy.md` §3.)
- **iVar** (2021, Genes) — direct prior art for "manage annotation versions → re-classify over
  time," but manual export/annotate/import round-trips, single-lab, no partitioning. VariantGrid
  does this internally with versioned partitions + automated diffs.
- **Exome reanalysis reviews / yield studies** (PMC8410709; Exomiser reinterpretation, npj Genom
  Med 2024; cancer variant reclassification PMC7469614) — quantify reanalysis payoff (≈10–15%
  extra yield; ~1/3 of reclassified variants change management). Justify the "genetic results have
  a long lifetime; automate reanalysis" argument.

### F. Scope decisions (what the VG paper does NOT centre)
- **Classification sharing / discordance → Shariant's turf** (Tudini 2022, AJHG). The VG paper
  says in one paragraph "VariantGrid is the technical platform beneath Shariant (ref)" and does not
  centre it.
- **HGVS resolution → cdot** (github.com/SACGF/cdot), to be **published first** and cited as a
  dependency. De-emphasise HGVS internals in the VG paper.

---

## 2a. Additional capability differentiators (from D. Lawrence, 2026-06)

- **Variant-centric AND sample-centric.** Users can search any variant/HGVS and, if absent,
  **create + annotate it on demand** — VariantGrid acts as a variant knowledge base
  (Variantopedia) as well as a sample-analysis tool. seqr and VarFish are organised around
  *loaded callsets*; an arbitrary variant not in a project cannot be looked up. Genuine
  capability difference.
- **Cross-build allele linking of annotations, tags, AND classifications.** Liftover links each
  variant to one `Allele` across GRCh37/GRCh38(/T2T) via the ClinGen Allele Registry, so a
  classification or tag created on one build is visible on the other. seqr/VarFish are typically
  single-build-per-project. Same backbone as the "one allele record" ingestion story (#5).
- **Routine clinical production + LIS integration APIs.** Deployed in accredited diagnostic use
  (~30k clinical samples) with laboratory-information-system integration — most academic tools
  never reach this. This, not raw record count, is the deployment credibility.

### Scale & deployment — frame honestly
- Dataset size is **mid-sized** (~5k research, ~30k clinical samples) vs seqr (23k exomes + 4k
  genomes / ~16,800 families) and OpenCGA (~100k genomes). **Do not lead with record count** —
  reviewers will just compare. Lead the scalability claim with an **architectural benchmark**
  (constant query time across annotation versions; SQL-array packing) and lead the maturity
  claim with **routine clinical deployment + LIS integration**.
- Setup is **non-trivial** (Postgres + Redis + RabbitMQ + Celery queues + VEP + hundreds of GB
  of annotation reference data). Planned **Docker image** is the mitigation to foreground. The
  honest differentiator is **on-premise on commodity hardware, no cloud/Kubernetes/big-data
  cluster** — not "easy to install."
- **Public vs self-hosted feature split** (affordability + venue, 2026-06): a public no-login
  instance can affordably offer **variant entry → on-demand VEP annotation → lookup →
  classification → share** (per-variant, cheap). **Node-based analysis on whole-genome/exome VCFs,
  ingestion, and the internal population DB stay self-hosted** (per-genome, expensive). This split
  satisfies NAR Web Server's no-login rule without giving away WGS compute (see
  `../journal_targets.md`). Caveat: the public slice overlaps free annotators (VarSome, Franklin,
  OpenCRAVAT, VEP-web) — its public hook must be **cross-build allele linking + shareable
  structured classifications**, not generic annotation.

## 3. Where VariantGrid is genuinely novel (defensible claims)

1. **Live interactive node-graph filtering.** Cartagenia Bench → Agilent **Alissa Interpret**
   pioneered node-based filter construction but runs **batch**; VarFish/seqr/Scout use form-based
   queries. **Alissa Interpret has been end-of-lifed** (Agilent announced discontinuation ~2023,
   retired end of 2024; migration path = QCI Interpret, which is NOT node-based). So VariantGrid
   is the **only actively-maintained node-graph variant analysis platform** — and the only
   open-source, live-editing one. *(Strongest UI novelty. The Alissa EOL strengthens this: state
   "only actively-maintained / only open-source live" rather than a bare "first.")*
   - **Source caveat:** Alissa EOL is from Agilent/QIAGEN vendor announcements + transition
     webinars (not a formal publication) — cite as a vendor EOL notice with accessed-date URL.
   - **Defensible claim (scope it to exploration):** node-graphs beat large filter forms *for
     open-ended research analysis where the right strategy is unknown*, because of (a) immediate
     per-edit feedback, (b) branching/merging of variant sets (Venn intersections, parallel
     strategies) that a linear form can't express, (c) inspectable intermediate results at every
     node, and (d) **the graph is its own provenance** — re-runnable, shareable, auditable; a
     filter form's state is ephemeral. Do NOT claim "better than filters, full stop" — that
     invites a counterexample. Keep the 2018-draft nuance: **interactive node-graph for
     exploratory research; batch/templates for fixed clinical pipelines** (VariantGrid does both,
     which is itself a strength).
   - **Evidence for the paper:** a worked case-study figure walking one real exploration (e.g.
     trio de-novo → gene list → Venn intersect → tag/classify) that explicitly shows a branch a
     linear filter cannot represent; plus the validated lineage of dataflow/visual programming
     (Make 1976; NUKE node compositing) brought *live* to variant analysis. Years of real use by
     dozens of scientists + reused templates = adoption evidence; a formal user study is likely
     out of scope.
   - **NAR demo angle:** make this feature publicly explorable via cached demo analyses on public
     cohorts (see `../journal_targets.md`) so reviewers experience it directly.
2. **Annotation-version partitioning for constant-time multi-version analysis + automated
   re-analysis via version diffs.** No comparator stores many full annotation versions with
   constant query cost and diffs them computationally. iVar and VarFish re-analyse but do not
   snapshot+diff at this granularity. *(Strongest architectural novelty.)*
3. **SQL-native multi-sample array packing** (Gemini extended into Postgres arrays so filtering
   stays in the query planner). *(Solid, concrete, benchmarkable contribution.)*
4. **End-to-end + multi-lab classification-sharing/discordance/ClinVar in one open platform**,
   evidenced at national scale by Shariant. *(The "whole lifecycle, deployed for real" story.)*

## 4. Where claims need softening / extra care

- "**First**" claims (live filters, partitioned variant DB): qualify to "first published /
  first open-source," and explicitly position vs **VarFish** (same stack) and **Alissa** (node UI).
- A **feature-comparison table** (the 2018 draft's "Supp Table 1") should now include VarFish,
  seqr, Scout, OpenCGA, GEMINI, VCF-Miner, Alissa, Franklin/VarSome — VarFish is the must-have row.
- Quantify the headline performance claims (constant-time annotation versioning; SQL-array
  filtering scaling) with a benchmark figure — reviewers at a methods journal will expect it.
- Treat **Shariant as a deployment/outcome of VariantGrid**, not a competitor, to avoid
  self-comparison confusion.

## 5. Inventory of retrieved literature

Full PDFs in `pdfs/` (→ text in `text/`): shariant_2022, vep_2016, gemini_2013, vcfminer_2016,
opencravat_2019, clinvar_2018, exome_reanalysis_review_2021, exomiser_reinterp_2024,
cancer_reclass_2020 (+ varsome/exomiser if recovered).
Web-extract text only (PDF auto-download blocked; see paywalled list): varfish_2020,
seqr_2022, ivar_2021, matchmaker_2015, decipher_2022, platform_compare_2025.
Paywalled / to retrieve manually: ACMG Richards 2015, ClinGen Rehm 2015, Sherloc 2017,
seqr final (Hum Mutat), reclassification clinical-impact studies — see `paywalled_articles.md`.
