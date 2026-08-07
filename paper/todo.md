# VariantGrid paper — TODO

Sequenced path to publication. Order reflects dependencies, not just priority.
Tags: **[you]** = you only (clinical env / submission / ethics) · **[claude]** = Claude can do/draft ·
**[both]**. See `strategy.md` for the why behind each.

---

## Now / in parallel (no blockers)

### 1. Upgrade SA Pathology to VG4  **[you]**
- [ ] Complete the VG4 upgrade (in a few months) → gives a **clinically-validated release**.
- [ ] After upgrade: migrate `vg3_sapath_prod` data up so VG4 data-mining scripts run on the full
      2020→present history.
- _Why first: unlocks both the "clinically validated, in accredited diagnostic production" claim
  and the data-mining (scripts target VG4 models)._

### 2. Pull out cdot and publish it  **[both]**  — companion paper, do BEFORE the VG paper
- [ ] Extract HGVS/transcript resolution into `cdot` (github.com/SACGF/cdot) as a standalone lib. **[you]**
- [ ] Write the cdot paper (short tool/Application note). **[claude can draft]**
- [ ] Submit cdot. **[you]**
- _Why early: removes HGVS from VG paper scope and gives a citable dependency. Cross-reference as a
  companion paper to pre-empt any salami-slicing concern._

### 3. Kiosk / cached public demo  **[both]**
- [ ] Stand up a **no-login cached demo**: pre-loaded analyses on public cohorts so anyone can drive
      the live node-graph (+ variant lookup/annotate/classify). **[you build / claude can help]**
- [ ] Optional: small **capped-VCF kiosk** (one ephemeral VCF). **[you]**
- _Why: unlocks **NAR Web Server** (no-login rule) and is a "try it now" credibility boost. Bounded
  scope — not a full SaaS. See `journal_targets.md`._

### 4. Start the data-mining governance path  **[you]** — LONG POLE, start ASAP
- [ ] Secure ethics/QA approval (HREC or QA-audit exemption) for aggregate SA Pathology usage mining.
- _Why now: the scripts/figures are fast; the approval is what gates the timeline._

---

## After VG4 migration + approval

### 5. Run the data-mining scripts  **[you]**
- [ ] Review `scripts/paper_data_mining.py`; confirm the `# VERIFY` fields against VG4 schema.
- [ ] Run read-only on a replica/backup → produce aggregate CSVs.
- [ ] Sanity-check 2-3 totals vs known VG3 values (migration preserved them).
- [ ] Bring de-identified CSVs back; Claude helps interpret + build figures. **[both]**
- _Target figures (highest value): ingestion + genotype:variant ratio · annotation-version diffs ·
  node-type usage + template-modification · classified-variant recurrence. See
  `data_collection_plan.md` / `data_mining_opportunities.md`._

---

## Paper-writing track  **[claude can draft, you finalise]**

### 6. Lock framing & venue
- [ ] Decide lead framing: software-architecture (→ NAR Web Server / Bioinformatics) vs clinical-
      usage-evidence (→ Genome Medicine). Depends on how strong the data-mining Results land.
- [ ] Confirm the NAR Web Server submission window if going that route.

### 7. Literature finishing
- [ ] Retrieve must-cite paywalled standards: **ACMG (Richards 2015)**, **ClinGen (Rehm 2015)** —
      see `literature/paywalled_articles.md`. **[you]**
- [ ] Verify the `?`/`P` comparator cells in `literature/feature_comparison.csv` vs current vendor/
      tool docs. **[both]**
- [ ] Condense the CSV to a **~15-row headline feature table** (STRONG cluster) + full CSV as supp. **[claude]**

### 8. Draft the manuscript
- [ ] Outline (intro → architecture → features → real-usage Results → discussion).
- [ ] Write architecture + figures (data model, annotation-version partitioning, node-graph).
- [ ] Configurability subsection + figure (one codebase → 4 deployments, feature-flag matrix).
- [ ] One-paragraph "substrate beneath Shariant (ref)"; cite cdot for HGVS; position Talos as
      complementary (don't out-yield).
- [ ] Honest "future work": phenotype ranking, Matchmaker Exchange.

### 9. Submit  **[you]**
- [ ] Internal review → submit → respond to reviewers.

---

## Parking lot (do NOT do for the paper unless a deployment needs it)
- Chasing feature parity (Exomiser-style phenotype ranking, Matchmaker integration) — low ROI;
  acknowledge as future work instead. _(See `strategy.md` §6.)_

## Critical path
VG4 upgrade ──▶ migrate data ──▶ run mining scripts ──▶ figures ──▶ submit
                  (ethics approval must be done before this step)
cdot paper + kiosk run in parallel and gate venue, not the data.
