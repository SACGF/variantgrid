# Data-mining opportunities — VariantGrid's unique empirical evidence

_Years of real clinical + research usage data is the asset NO competitor's paper can replicate.
A software paper backed by real-world usage analytics is rare and persuasive. This is likely the
**empirical heart / Results section** of the VariantGrid paper. Captured from D. Lawrence's list
plus extensions._

## Why this matters
Most tool papers prove value with one curated cohort or a yield number. VariantGrid has been in
**routine accredited clinical diagnostic use (SA Pathology) and research use for years** — so it
can show **how variant interpretation is actually done in practice**, longitudinally. That is
evidence, not features, and it is unique to a long-lived deployed platform.

## Minable questions (grouped)

### A. Growth & re-analysis over time (ties to the annotation-versioning story)
- VCFs / samples / variants / genotypes ingested over time (the genotype:variant ratio curve —
  the 2018 draft's "668 exomes → 2.4B genotypes vs 46M variants, ~50x"; refresh with current
  numbers).
- Novel-variants-per-VCF decay curve as the database grows (justifies managed shared annotation).
- How often annotation versions changed, and **what changed between versions** (use the
  version-diff transition matrices) — e.g. ClinVar significance churn, gnomAD AF shifts.
- Cases where a stored classification's underlying evidence changed after the fact (the
  reanalysis payoff, VariantGrid-native — complements Talos's yield numbers).

### B. How people actually do analysis (unique to the node-graph)
- Distribution of node types used; typical analysis graph size/shape/depth.
- Which nodes users add/remove/reorder during a session (the "live editing" value, measurable).
- What users **modified from the source/template** — how often templates are tweaked vs run as-is
  (supports "interactive exploration vs fixed clinical pipeline" thesis with real data).
- Time-to-first-candidate / iterations per analysis.
- Most-used columns out of the 100+ available.

### C. Tagging & curation behaviour
- What variants get tagged, with which tags, and the funnel from tagged → classified → reported.
- How often a tagged/curated variant recurs across later samples (internal-population-DB value).
- Lag between first-seen and first-classified for a variant.

### D. Where classified variants turn up (internal matching value)
- How often a newly sequenced sample carries an already-classified variant (auto-flag value).
- Variants classified in one lab/sample later seen in another (internal matchmaking).
- Founder / recurrent variants surfaced by internal frequency that gnomAD under-represents
  (artifact-vs-real story).

### E. Search behaviour
- What users search for (HGVS vs gene vs coordinate vs rsID vs phenotype) — note: HGVS resolution
  itself is moving to **cdot**, but search *usage patterns* stay a VariantGrid result.
- Search → create-on-demand-annotation frequency (how often a searched variant was novel).

### F. Multi-deployment / configurability evidence
- Contrast usage profiles across the 4 deployments (research vs clinical diagnostic vs Shariant
  vs RUNX1db) from the SAME codebase — shows the platform flexes to very different workflows.

## Cautions
- **Ethics/consent/governance:** clinical-usage mining needs appropriate approvals; aggregate,
  de-identified analytics only. Confirm before any figures leave the building.
- Keep patient-level data out; report rates, distributions, curves — not records.
- Some of this is its own potential paper; for the platform paper, pick 3-4 figures that directly
  back the architecture claims (A and B are highest-value).

## Highest-value figures for the paper (recommendation)
1. **Genotype:variant ratio + novel-variant decay** over ingestion time → backs SQL-packing +
   managed-annotation scaling claims.
2. **Annotation-version change matrix** (what flipped between two real versions) → backs the
   versioned-annotation + reanalysis claim natively.
3. **Node-type usage + template-modification rates** → backs the live node-graph exploration claim.
4. **Already-classified-variant recurrence rate** in new samples → backs internal population DB /
   matching claim.
