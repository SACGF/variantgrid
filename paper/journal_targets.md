# VariantGrid paper — candidate journals

_Rationale and ranked targets, given the paper is an **open-source clinical-genomics software/
platform** paper whose novelty is architectural (annotation-version partitioning, live node-graph
filtering, SQL-native multi-sample packing) plus a real-world national deployment (Shariant)._

## Audience (drives abstract framing more than journal choice)

Two-tier:
- **Champions / deciders — lab heads, clinical lab directors, pathologists.** Decide to adopt.
  Care about: diagnostic yield + automated re-analysis, inter-lab discordance/quality (Shariant),
  **data sovereignty** (self-hosted, patient data never leaves the institution), free-for-research,
  ClinVar participation, audit trail.
- **Implementers — the lab's own bioinformatician / sysadmin.** Deploy and maintain. Care about:
  familiar stack (Django/PostgreSQL), install + maintenance burden, scalability, APIs/hooks.

Goal stated by D. Lawrence: *get lab heads to have their bioinformaticians deploy it.* So the
paper must make the lab head **want** it (abstract + intro lead with clinical/quality/sovereignty
value) while giving the bioinformatician confidence it is deployable.

**Honest deployment framing (do NOT overclaim "trivial / one person").** VariantGrid is a
multi-service stack — PostgreSQL, Redis, RabbitMQ, multiple Celery worker queues, VEP + its
large annotation reference data (hundreds of GB). It is **non-trivial to set up**; a Docker
image is planned and should be foregrounded *because* of that, not to imply it is easy. The
true, defensible adoption advantage is **self-hostable on-premise on commodity hardware — no
managed cloud, no Kubernetes, no Hadoop/Spark/Solr cluster** — vs seqr (GCP/Kubernetes) and
OpenCGA (big-data stack). (VarFish is the fair comparator: similarly heavy stack, also ships
container/Snakemake deployment.) State "self-hostable without a cloud or big-data team," keep
the data-sovereignty argument, and let the Docker image carry the ease-of-deployment message.
NAR Web Server / Bioinformatics are read by exactly these implementers, so the venue is right
**provided the abstract sells lab-head value**.

## What kind of paper this is

A **software / application** paper with (a) a genuine methods contribution (the database
architecture and re-analysis design) and (b) substantial deployment evidence. That profile fits
either a **bioinformatics-software** venue or an **application-notes** track. The comparators
landed in a consistent set of journals, which is the strongest signal for where this belongs:

- VarFish, ANNOVAR → **Nucleic Acids Research (Web Server / Database issue)**
- VEP → **Genome Biology**
- GEMINI → **PLOS Computational Biology**
- seqr, DECIPHER, Matchmaker → **Human Mutation** (now **Human Genetics and Genomics Advances**)
- Shariant → **American Journal of Human Genetics**

## Deciding constraint: public no-login access (raised by D. Lawrence, 2026-06)

**NAR Web Server requires the server be freely usable with NO login**, and must stay available
≥2 years. VariantGrid is primarily self-hosted and cannot give the world free whole-genome VCF
upload + node-based analysis.

**Resolution (updated 2026-06): split the product into two layers.** D. Lawrence confirms a
public no-login instance CAN affordably offer **variant entry → on-demand VEP annotation →
variant lookup → classification (+ browse/share)** — cheap because it is per-variant, not
per-genome. The expensive node-graph analysis on whole-genome/exome VCFs + continuous ingestion
+ internal population DB stays **self-hosted/downloadable**.

This split maps onto the classic NAR Web Server model (permanent public server for the
interactive part + downloadable for heavy/local use). So **no throwaway "kiosk" is needed** — the
locked-down public instance is a real, permanent, useful resource, and **NAR Web Server is viable
again as a primary option** without giving away WGS compute. (VarFish precedent still applies; we
simply don't need its kiosk hack.)

**Caveat that drives framing:** the publicly-exposed slice (annotate one variant + classify) is
the LEAST differentiated part of VariantGrid — it overlaps with VarSome, Franklin, OpenCRAVAT and
Ensembl VEP-web, which annotate variants free. For NAR the public instance's hook must be what
those tools do NOT do:
- **Cross-build allele linking** — a classification/tag made on GRCh37 is visible on GRCh38 and
  vice-versa (ClinGen Allele Registry backbone).
- **Structured, shareable classifications** with the evidence-key model + Shariant sharing story —
  a classification *knowledge base*, not a one-shot annotator.
Lead the public framing on those two and it carries an NAR paper; lead it on "we annotate
variants" and it does not.

**Stronger option (D. Lawrence, 2026-06): cached/capped public node-graph demo.** A no-login
instance can ship **pre-loaded demo analyses on public cohorts** (existing CCB research exomes on
variantgrid.com, or a public 1000G/platinum set) so reviewers interactively **click nodes,
re-filter, see counts, drill into grids** — experiencing the headline live node-graph feature at
near-zero marginal compute (nothing uploaded/freshly annotated). A small **capped-VCF kiosk** adds
one fresh single-case run. Tiering: (1) no-login cached exploration + lookup/annotate/classify;
(2) no-login capped kiosk (one small VCF); (3) self-hosted unlimited ingestion + WGS + population DB.

**This resolves the earlier "public slice is least novel" caveat** — the public instance now
demonstrates the MOST novel capability (live node-graph exploration), not just generic annotation.
**NAR Web Server therefore moves from "viable" to the clear primary**, and the framing tension
below mostly dissolves (we no longer must lead on annotation to satisfy the no-login rule).

## Ranked targets

_#1 vs #2 is now a FRAMING choice (see "Deciding constraint" above), not a feasibility gate —
both are viable. Pick by which story leads._

### 1a. Nucleic Acids Research — Web Server issue  (if public interpretation layer leads)
- **Why:** natural home for a maintained genomics web server; VarFish (closest comparator) and
  ANNOVAR are here; high visibility; open access; ideal reviewer community.
- **Fit:** viable now that the affordable public layer (annotate → lookup → classify → share,
  cross-build) satisfies the no-login rule without giving away node-based WGS compute. Hook the
  public instance on **cross-build allele linking + shareable classifications**, not generic
  annotation (which VarSome/VEP-web already do free).
- **Watch-outs:** annual submission window (deadlines early in the year); public instance must be
  live + useful + maintained ≥2 years (variantgrid.com already runs, so OK); reviewers must see a
  hook beyond "annotate a variant."

### 1b. Bioinformatics (Oxford) — full Application paper  (if self-hosted architecture leads)
- **Why:** standard venue for bioinformatics software; broad implementer readership; cyvcf2, vt,
  VarSome appeared here. A full paper accommodates the architecture + benchmarks + deployment
  narrative.
- **Fit:** strong for the methods/engineering angle; **no mandatory no-login public server** —
  freely available source + Docker image + a demo satisfies availability. Best home for the
  node-graph / partitioning / population-DB story.
- **Watch-outs:** slightly less clinical-genetics framing than AJHG/HGGA; make the availability
  statement (source + Docker) prominent.

### 3. GigaScience  (or GigaByte)
- **Why:** explicitly rewards reproducible, large-scale, open data/software with living
  resources; would welcome the scalability benchmarks and open deployment. Good if we want to
  emphasise the data-engineering contribution and provide reproducible benchmark artifacts.
- **Fit:** very good for the "scalable open platform" framing; lower clinical readership.

### 4. Genome Medicine
- **Why:** if we lead with **clinical impact + Shariant deployment** (re-analysis yield,
  inter-lab discordance resolution at national scale) rather than pure software. Clinical
  genomics audience.
- **Fit:** good if the framing is "platform enabling clinical outcomes"; needs outcome data to
  carry it.

### 5. Human Genetics and Genomics Advances (HGGA, open-access sister of AJHG)
- **Why:** seqr, DECIPHER, Matchmaker and the Shariant work live in this AJHG/Human-Mutation
  family. Right audience (clinical/statistical genetics, data sharing). Open access.
- **Fit:** good, especially to keep continuity with the Shariant (AJHG 2022) paper and the
  data-sharing framing.

### Also reasonable
- **BMC Bioinformatics** (solid, lower-impact software home).
- **JAMIA / J Biomedical Informatics** (if the LIS-integration / informatics-system angle leads).
- **F1000Research / Wellcome Open** (fast, versioned, open — fallback or companion).

## Recommendation

Both **NAR Web Server** and **Bioinformatics** are now viable (the public/self-hosted layer split
removes the access blocker). Choose by framing:
- **NAR Web Server** if we lead with the publicly-accessible interpretation/classification layer —
  hooked on cross-build allele linking + shareable classifications. Higher visibility, ideal
  reviewers, VarFish precedent. Requires maintaining the public instance ≥2 years.
- **Bioinformatics (full paper)** if we lead with the self-hosted architecture (node-graph,
  annotation-version partitioning, internal population DB). No public-server obligation.

My lean: **NAR Web Server (primary).** With a cached/capped public instance that lets reviewers
interactively explore the live node-graph, it is the higher-visibility venue in the exact
community AND publicly demonstrates the headline feature — no WGS compute given away. The main
remaining commitment is maintaining the public instance ≥2 years (variantgrid.com already runs).
**Bioinformatics** stays the safe, equally-respectable fallback if we'd rather not run a public
instance. For a clinical-outcomes lead use **Genome Medicine / HGGA**; for a reproducible-benchmark
lead, **GigaScience**.

Two framing decisions drive the choice and are worth settling before drafting:
1. **Software-architecture paper** (lead with the DB design + node-graph + benchmarks) → NAR /
   Bioinformatics / GigaScience.
2. **Clinical-platform paper** (lead with Shariant deployment + re-analysis yield) → Genome
   Medicine / HGGA.
The architecture framing is the more defensible novelty and has the cleaner comparator set, so
NAR Web Server is the lead recommendation.
