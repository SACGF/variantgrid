# VariantGrid paper — strategy guide

_Synthesis of the positioning discussion (2026-06). Answers: what angle, what's in/out, and where
to spend effort (feature completeness vs kiosk vs data-mining). Companion files: `journal_targets.md`,
`data_mining_opportunities.md`, `literature/comparison.md`, `literature/feature_comparison.csv`._

---

## 1. The angle (thesis)

> **VariantGrid is a configurable, clinically-validated, self-hostable platform whose data model
> enables continuous ingestion, versioned annotation, and automated re-analysis of genomic
> variants — with interactive node-graph analysis — proven by years of real-world clinical and
> research deployment.**

Lead with the **platform + data model + real usage evidence**. Treat sharing/discordance and HGVS
as substrate you provide, not the focus (they have their own homes — see scope). The differentiated,
hard-to-replicate assets are the architecture, the configurability, the clinical validation, and —
most of all — the **real usage data** (no competitor's paper has this).

Don't frame it as "the biggest" or "the most features." Frame it as **"maturely engineered,
deployed for real, and architected for the long life of genetic data."**

---

## 2. Scope — what's IN and what's OUT

| Topic | Decision | Why |
|-------|----------|-----|
| Continuous ingestion → single cross-build allele record | **IN (lead)** | Core novelty; enables internal population DB + matching |
| Versioned annotation (partition-per-version) + version-diff re-analysis | **IN (lead)** | Strongest architectural novelty; nobody else does it |
| Live interactive node-graph analysis | **IN (lead)** | Only actively-maintained node-graph tool (Alissa EOL) |
| Internal population DB + sample-level "who carries this" | **IN** | Clean differentiator vs seqr; pairs with allele model |
| **Configurability** — one codebase → 4 deployments | **IN (elevate)** | New, distinctive; see §4 |
| **Clinical validation** (VG4 at SA Pathology) | **IN (elevate)** | Credibility most academic tools never reach; see §5 |
| **Real-world usage data-mining** | **IN (the empirical Results)** | Unique evidence; see `data_mining_opportunities.md` |
| Classification sharing / discordance | **OUT → brief** | **Shariant's turf** (Tudini 2022, AJHG). Say "VariantGrid is the technical platform beneath Shariant (ref)" in one paragraph; do not centre it |
| HGVS resolution | **OUT → separate** | Moving to **cdot** (github.com/SACGF/cdot); **publish cdot FIRST**, then cite it as a dependency. Removes scope + gives a citable building block |
| ACMG classification mechanics | **IN (light)** | Mention as a capability; the deep evidence-key/discordance detail belongs with Shariant |
| Phenotype ranking / Matchmaker | **OUT (acknowledge)** | VG is behind here; name it as out-of-scope / future work, don't overclaim |

---

## 3. Talos and the re-analysis positioning (important nuance)

**Talos** (Centre for Population Genomics, medRxiv 2025; `literature/talos_2025.txt`) is the current
state-of-the-art **automated re-analysis** tool — Australian, open-source, **batch pipeline**
(Nextflow/Hail/Docker, cohort-scale, monthly ClinVar/PanelApp cycles, JSON output), 5.2% additional
yield on 4,735 undiagnosed patients. It compares itself to Exomiser.

**Do NOT try to out-yield Talos.** It is purpose-built for automated reanalysis yield and is recent
and strong. VariantGrid's reanalysis story is **different and complementary**:
- Talos = a standalone automated **flagging pipeline** (no UI, no persistent queryable database, no
  interactive analysis).
- VariantGrid = a **platform** where reanalysis is *integrated* — annotation-version diffing inside a
  queryable database that already holds the samples, classifications, tags, and analyses; plus
  interactive node-graph exploration Talos doesn't offer.
- They could **interoperate** (Talos as an engine; VariantGrid as the platform/substrate).

Frame VG as **infrastructure for reanalysis** (versioned annotation + diffs + storage + UI), and
cite Talos + iVar + Exomiser-reinterpretation + the reanalysis-yield reviews as establishing that
reanalysis *matters*. Position, don't compete.

---

## 4. The configurability differentiator (elevate this)

One codebase, feature-flagged + lightly reskinned, runs as **four very different production
systems**: Shariant (national classification-sharing network), RUNX1db (disease-specific variant
database), a public research server (variantgrid.com), and an **accredited clinical diagnostic
server** (SA Pathology). Source-backed: hostname-based split settings + feature flags
(`UPLOAD_ENABLED`, `DISCORDANCE_ENABLED`, `CLINVAR_EXPORT_ENABLED`, `SHARIANT_SYNC_ENABLED`, …).

Why it's worth elevating:
- It's genuinely unusual — most tools are one deployment shape.
- It speaks directly to the **lab-head audience**: "you can run it as *your* thing" (your branding,
  your enabled feature set, your data).
- It's evidence of engineering maturity and of a clean architecture.

Add a row to the feature table and a short subsection + a figure (the 4 deployments from one
codebase, with the feature-flag matrix).

---

## 5. Timing & sequencing

1. **Publish cdot first** (HGVS → github.com/SACGF/cdot). Clears HGVS out of the VG paper scope and
   creates a citable dependency. Lower-effort, standalone, de-risks the bigger paper.
2. **Time the VG paper to the VG4 / SA Pathology upgrade** so the paper can state *"clinically
   validated; in routine accredited diagnostic production (VG4)."* That sentence is worth more than
   any feature row.
3. **Instrument + extract the usage-mining figures** (needs governance sign-off; start the ethics/
   aggregation path early — it's the long pole).
4. **Stand up the cached public demo** (see §6) in parallel — cheap given high dev velocity.

---

## 6. THE STRATEGIC QUESTION: where to spend your time

You asked: chase feature completeness, or build a kiosk / public login-free version?

**Answer: neither as the priority — spend it on (a) the data-mining results and (b) a cheap cached
public demo. Skip feature-completeness chasing.**

### Do NOT chase feature completeness
- Reviewers do not reward parity ("we also have X"). The feature table is **defensive** (proves you
  miss no basics), not the story.
- The two gaps (Exomiser-style phenotype ranking, Matchmaker) are real but **low ROI** — acknowledge
  as future work. Only build them if a *deployment* needs them, not for the paper.
- Your STRONG cluster is already large and differentiated. Adding more features dilutes effort.

### DO spend time on the data-mining results (highest ROI)
- It's **evidence no competitor can produce**, it plays to your high dev velocity + Claude Code, and
  you already have the data. It converts a "here's our software" paper into "here's how genomic
  interpretation actually happens at scale, over years."
- Start the governance/ethics + aggregation path now (the slow part). Target 3-4 figures (see
  `data_mining_opportunities.md` §"highest-value figures").

### DO build a cached public demo (high ROI, bounded effort)
- A no-login instance with **pre-loaded demo analyses on public cohorts** lets reviewers/readers
  interactively explore the live node-graph (your headline feature) at near-zero compute — plus
  variant lookup/annotate/classify. Optional small capped-VCF kiosk on top.
- This is what **unlocks NAR Web Server** (no-login requirement) AND is a "try it now" credibility
  boost. Bounded scope — don't build a full multi-tenant SaaS; the affordable interpretation layer +
  cached demos is enough.

### Kiosk vs data-mining, if you must choose one
- **Data-mining wins for paper strength** (it's the differentiated evidence).
- **Kiosk wins for venue + adoption** (unlocks NAR, lets people try it).
- Given your dev velocity the cached demo is cheap, so realistically **do both** — data-mining is the
  content, the demo is the shop window.

---

## 7. Venue implication (see `journal_targets.md` for detail)
- The **data-mining + clinical-validation** substance can push toward a higher-impact clinical-
  informatics framing (**Genome Medicine / HGGA**) OR strengthen an **NAR Web Server** submission
  (if the cached demo is built) / **Bioinformatics** (if not).
- My lean: build the cached demo → **NAR Web Server** primary (highest visibility, ideal reviewers,
  headline feature publicly demonstrable). If the public-instance commitment is unattractive,
  **Bioinformatics**. If the usage-data Results become the star, consider **Genome Medicine**.

---

## 8. Recommended next actions (prioritised)
1. **Start the usage-data-mining governance + extraction** (long pole; unique value). Pick the 3-4
   figures from `data_mining_opportunities.md`.
2. **Publish cdot** (clears HGVS scope; citable dependency).
3. **Stand up cached public demo analyses** (unlocks NAR; "try it now").
4. **Time submission to the VG4 / SA Pathology clinical-validated release.**
5. **Write the ~15-row headline feature table** (from the STRONG cluster) + keep full CSV as supp.
6. **Verify the `?`/`P` comparator cells** in the CSV against current vendor/tool docs.
7. Retrieve the must-cite paywalled standards (ACMG 2015, ClinGen 2015) — see
   `literature/paywalled_articles.md`.

## One-paragraph summary
Lead with the **platform and its data model** (continuous ingestion → one allele record, versioned
annotation with diff-driven re-analysis, internal population DB), the **live node-graph** (now the
only actively-maintained one), the **configurability** (one codebase → 4 deployments), and the
**clinical validation** (VG4 at SA Pathology) — and make the **empirical Results a data-mining of
real usage** that no competitor can match. Push sharing/discordance to a one-paragraph "we are the
substrate beneath Shariant" and move HGVS to a separately-published cdot. Don't chase feature parity;
spend your velocity on the usage-data figures and a cheap cached public demo. That demo unlocks NAR
Web Server; without it, Bioinformatics; if the usage data steals the show, Genome Medicine.
