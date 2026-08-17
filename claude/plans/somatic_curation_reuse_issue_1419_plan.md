# Reusing prior curation for somatic reporting

Design for [#1419](https://github.com/SACGF/variantgrid/issues/1419) (separate out gene / disease
curation), [sapath#246](https://github.com/SACGF/variantgrid_sapath/issues/246) (SomaticReportable →
easy classifications), and the entry point [#444](https://github.com/SACGF/variantgrid/issues/444)
will consume. This is Phase 8 of [`tso500_overall_plan.md`](tso500_overall_plan.md).

The driving observation: in somatic work the same variants recur constantly, and the lab has usually
curated them before. Copy consensus already reuses prior work per allele. What is missing is a way to
*see* what has been curated before you start, a way to get through a report's worth of variants without
re-navigating for each one, and a way to reuse the gene-level content that is identical across every
variant in that gene.

---

## What already exists

Most of the machinery is here; the delta is smaller than the issue titles suggest.

- **Allele-level copy consensus.** `ClassificationConsensus.all_consensus_candidates()`
  (`classification/models/classification.py:2664`) returns *Latest Germline* and *Latest Somatic* for an
  allele — published, non-withdrawn, `exclude_external_labs=True` — defaulting the radio to whichever
  matches the user's `allele_origin_focus`. `consensus_patch` (`:2685`) copies `value` and `note` for
  every ekey with `copy_consensus=True`, and forces `allele_origin` from the source's bucket rather than
  copying it.

- **Somatic vs germline is data, not code.** The `allele_origin` option `somatic` carries
  `namespaces: ["somatic"]` and `assertion_method` option `amp` carries `["amp"]`;
  `_evidence_key_overrides_from_evidence_fields` (`classification/models/classification.py:1965`) turns
  those into the per-record namespace set that switches the `amp:`/`horak:`/`somatic:` keys on. So
  "launch as AMP/somatic" means seeding those two keys at creation.

- **Tag → classification.** `CreateClassificationForVariantTagView` (`analysis/views/views.py:1207`)
  already narrows the sample dropdown to the tag's node samples and posts to
  `create_classification_for_analysis`, which links the new record to the analysis via
  `AnalysisClassification` (`:1275`). The only reason the path is germline-flavoured is
  `formatVariantTagFirstColumn` (`variantgrid/sitestatic/static/js/grid.js:541`), which hardcodes the
  New Classification button to `RequiresClassification`.

- **Gene lookup.** `ResolvedVariantInfo.gene_symbol` is a real FK
  (`classification/models/classification_variant_info_models.py:143`) and
  `classification_gene_symbol_filter` exists, so "prior classifications in this gene" is one query with
  no new denormalisation.

The gap on the specimen side: `VariantTag` (`analysis/models/models_variant_tag.py:42-51`) has no
sample. Tag → specimen goes tag → `node.get_samples()` → `Sample.extraction` → `Specimen`.

---

## Part A — the copy consensus audit

142 of 330 evidence keys have `copy_consensus=True`. The flag was set for germline ACMG work and has
never been reviewed against somatic. Four findings, in descending severity.

The audit is split across two issues, both raised, and both worth checking the status of before
starting the rest of this plan:

- **[#1713](https://github.com/SACGF/variantgrid/issues/1713)** — A1, A2 and A3's namespace filter. Wrong
  today, needs no new vocabulary, wanted ASAP.
- **[#1714](https://github.com/SACGF/variantgrid/issues/1714)** — the `copy_scope` / `copy_allele_origin`
  fields, covering A3's twenty germline-only keys and Part B's gene scope. Goes through Shariant triage,
  since it changes a field labs configure.

Land #1713 first: it is independent, and its `copy_consensus=False` settings migrate straight to
`copy_scope = NONE` when #1714 follows.

### A1 — Test-level facts about *this* tumour are marked copyable

| Key | Category | Why it is wrong |
|---|---|---|
| `somatic:tmb_value`, `somatic:tmb_status` | HT | TMB of the specimen being reported |
| `somatic:msi_value`, `somatic:msi_status` | HT | MSI of the specimen being reported |
| `somatic:hrd_status` | HT | HRD of the specimen being reported |
| `testing_context` | HT | Per-test |

These describe the patient's tumour, not the allele. Copying them from a prior classification imports
another patient's measurements. `somatic:tumor_cellularity` is already `False`, so the set is internally
inconsistent as well as wrong.

Phase 4's `SpecimenMeasure` (#1559) is the correct source for the five `somatic:` ones — they should
arrive through `evidence_from_sample_and_patient.py`, not through a copy.

**Set `copy_scope = NONE`.**

### A2 — Report and reviewer free text is marked copyable

`somatic:summary_interpretation` is described in its own migration
(`classification/migrations/0164_somatic_hrd_msi_tmb_ekeys.py`) as a free-text summary of all
classifications in a report, `max_share_level='lab'`, flagged there as very likely to contain
patient-identifiable information. Copy consensus pulls another patient's report narrative into this
record.

`review_comment` (SO) has the same shape — a reviewer's comment about a different record.

**Set `copy_scope = NONE` for both.**

### A3 — Cross-bucket leakage

`consensus_patch` uses `EvidenceKeyMap.cached()`, the unconfigured global map, so it copies all 142 keys
regardless of the target record's namespaces. `exclude_namespace` (`evidence_key.py:567`) only
suppresses mandatory-field validation (`classification.py:1241`) — it does not stop a patch writing the
key. So copying a germline record into a somatic one writes 28 `acmg:*` criteria into evidence that the
form will never show, and the reverse writes `amp:`/`horak:`/`somatic:` into a germline record.

Namespaced breakdown of the 142: 87 un-namespaced, 28 `acmg`, 17 `horak`, 6 `somatic`, 4 `amp`.

**Fix: filter the consensus patch by the target classification's `evidence_key_overrides.namespaces`.**
This is 55 of the leaky keys handled with no per-key data entry, and it stays scoped to the UI-initiated
copy — API imports keep writing out-of-namespace keys deliberately, which is how records shared from
labs with a different config survive.

That leaves germline concepts that have no namespace to filter on, all currently copyable:

```
condition_incidence  mode_of_inheritance  gene_penetrance  variant_penetrance  proband_count
segregation  segregation_meioses  segregation_affectedcarriers  segregation_unaffectedcarriers
segregation_affectednoncarriers  segregation_bayes  segregation_lod  s_other  s_summary
denovo_points  d_other  d_summary  a_other  a_summary  match_maker_exchange
```

**Fix: an explicit `copy_allele_origin` on `EvidenceKey` (ANY / GERMLINE / SOMATIC), set to GERMLINE for
the twenty above.** Twenty rows of data entry rather than a rule, because there is nothing in the schema
that distinguishes them.

### A4 — Redundant but harmless, leave alone

The gnomAD count keys and about a dozen predictor scores (`alphamissense`, `bayesdel`, `vest`,
`varity_r`, `mpc`, `primateai`, `metarnn`, `clinpred`, `mutpred2`, `aloft`) are copyable while their
siblings (`revel`, `cadd`, `sift`, `polyphen`) are not. Inconsistent, but autopopulate runs first and
wins: `create_classification_object` applies the consensus patch with `leave_existing_values=True`
(`classification/views/views.py:428-440`) and `AutopopulateView` guards on `used_keys` (`:349-356`), so
a copied value only ever fills a gap annotation left. Not worth churning.

### The `condition` key

`condition` is copyable today. For somatic it holds the tumour type, which is a per-patient fact and the
single axis on which reuse is most likely to go clinically wrong. It stays copyable, but the wizard
surfaces it as an explicit human decision rather than a silent pre-fill (Part C, decision 3).

### Model change — #1714

Replace the `copy_consensus` boolean with `copy_scope` (`NONE` / `ALLELE` / `GENE`) and add
`copy_allele_origin` (`ANY` / `GERMLINE` / `SOMATIC`). Migration: `True` → `ALLELE`, `False` → `NONE`,
then apply A1/A2's `NONE`, A3's twenty `GERMLINE`, and Part B's `GENE` set. Both fields belong in the
admin fieldset beside `max_share_level` (`classification/admin/classification_admin.py:640`), and
`legacy_somatic.py:287` follows.

---

## Part B — gene-level reuse, with the human picking

The gene-level content of a somatic classification — what the gene does, its role in cancer, the
gene-level literature — is identical across every variant in that gene, and is exactly what
`copy_consensus` cannot reach today because it is keyed on the allele.

### The candidates are shown, never auto-picked

AMP tiering and therapy content are gene **and tumour type**, not gene. Copying gene-level content from
a colorectal case onto a melanoma case is a clinical error rather than staleness, and the phenotype data
available to match on is not good enough to automate the judgement. So the gene-level source is always a
human choice.

The query: `ClassificationModification.latest_for_user(user, published=True)` filtered through
`classification_gene_symbol_filter(gene_symbol)` and the current allele-origin bucket, most recent
first, excluding whatever is already offered as an allele-level candidate, capped at ten. Each row shows
condition, clinical significance / tier, lab and curated date — enough to judge tumour-type relevance at
a glance. "None" is a first-class choice and the default.

### Which keys travel at gene scope

The `H` (Gene) category, less `condition` itself:

```
condition_incidence  disease_onset  essential_gene_crispr  essential_gene_crispr2
essential_gene_gene_trap  gene_constraint  gene_damage_index_score  gene_disease_validity
gene_indispensability_score  gene_penetrance  ghis  gnomad_oe_lof  gnomad_pli  gnomad_pnull
gnomad_prec  h_summary  hipred_score  mechanism_of_disease  mode_of_inheritance  phi  prec
variant_penetrance
```

plus `pubmed_gene_search_count` (L). Several of these carry `copy_allele_origin = GERMLINE` from A3,
which is correct: they are gene-level *and* germline-only, and both filters apply.

`literature` stays at `ALLELE`. It currently mixes gene-level and variant-level content, and splitting
it is [variantgrid_private#1102](https://github.com/SACGF/variantgrid_private/issues/1102) — this plan
waits for that rather than copying variant literature across a gene. `search_terms` stays at `ALLELE`
for the same reason: it contains variant terms.

### Merge order

Autopopulate wins over allele-level, which wins over gene-level. Allele beating gene is the specific
beating the general; autopopulate beating both is existing behaviour and stays.

### Other labs

Show them, do not copy from them. `all_consensus_candidates` passes `exclude_external_labs=True` today,
which hides what other labs — including Shariant — have curated on the same allele, and that is worth
seeing when deciding how to curate. But an external record was curated under another lab's config,
assertion method and namespaces, so copying its evidence into a local record imports assumptions that
were never reviewed here. External candidates therefore render in the overview as read-only context,
without a copy control.

### Not doing: a first-class gene/disease object

#1419's schema-change option — gene/disease curation as its own model with its own review cycle — is the
right eventual answer, and is the only thing that can express "reviewed on this date, due for review".
It needs the disease axis settled first (condition matching, ontology hierarchy, whether a gene/disease
record is per-lab), which is a much larger piece of work. The copy route above delivers the reuse now
using machinery that already works, and does not foreclose the object later: `copy_scope = GENE` is
exactly the set of keys such an object would own.

---

## Part C — the wizard

A report's worth of somatic variants is a batch, and the current flow makes it N independent
navigations. But each variant still needs real decisions, so batching the *decisions* would trade
navigation cost for a screen nobody can reason about. The wizard batches the triage and then serialises
the decisions, one variant at a time.

### Stage 1 — triage the list

One row per `SomaticReportable` tag in scope. Each row shows the variant and its gene, and a brief
overview of what already exists for that allele: how many prior classifications, the latest one's
clinical significance / tier, its condition, lab and date — external labs included, marked as such.

The scientist drops the rows they are not going to report. Dropping removes the `SomaticReportable` tag,
which is the same gesture as untagging in the analysis — a plain delete, with the existing post-delete
signal (`analysis/signals/signal_handlers.py:29`) updating node counts and nothing else recorded. Then:
*Classify N variants*.

Soft-deleting `VariantTag` was considered and rejected. There are 32 query sites across 15 modules
(grids, node counts, `get_for_build`, liftover tasks, the analysis JS tag dict), every one of which would
have to learn to exclude deleted rows, and `get_or_create` in `set_variant_tag`
(`analysis/views/views_json.py:227`) would need to undelete rather than duplicate. That is a lot of blast
radius on a high-frequency mundane gesture.

Classifying does **not** remove the tag — sapath#246 is explicit that the somatic lab uses
`SomaticReportable` to track somatic variants, so the tag survives the classification that came from it.

### Stage 2 — one variant at a time

For each kept variant in turn, a single screen carrying only the decisions:

1. **Which prior classification for this allele to copy from** — radio over the candidates, each showing
   condition, clinical significance / tier, lab, date. Plus "none".
2. **Which prior classification for this gene to copy from** — same shape, only shown when the gene has
   candidates, defaulting to none (Part B).
3. **The condition for this record** — pre-filled from whichever record was picked in 1 or 2, always
   visible and always editable.

Lab, sample and transcript resolve silently: lab from user settings, sample from the tag's node
(existing `_get_sample_form` behaviour, `analysis/views/views.py:1223`), transcript from the canonical
selection. An "advanced" disclosure exposes all three for the cases where the default is wrong. Create,
and the wizard advances to the next variant.

`allele_origin=somatic` and `assertion_method=amp` are seeded at creation, which is what makes the
record AMP/somatic (Part A's namespace mechanism). Copying from a somatic record already sets
`allele_origin` via `consensus_patch`; seeding it explicitly covers the "none" case.

### State

There is none to store. The tags are the worklist, and a tag is "done" when a classification for that
variant is linked to the analysis through `AnalysisClassification`. That makes the wizard resumable by
construction, survives a browser crash, and lets two people work the same list without a lock. Progress
("4 of 11") is a count, not a record.

---

## Part D — where it launches from

**A `TAG_SOMATIC_REPORTABLE` setting**, mirroring `TAG_REQUIRES_CLASSIFICATION`
(`variantgrid/settings/components/default_settings.py:720`), defaulting to `None` so the feature is off
except where a deployment names its tag. `formatVariantTagFirstColumn` stops hardcoding
`RequiresClassification` and drives off the configured tags instead.

**From the analysis**, beside the existing tags button — "Classify somatic reportable (N)" — scoped to
that analysis's tags.

**From the specimen, extraction and patient pages**, scoped to the tags on variants in analyses
containing a sample of that specimen's extractions. Two things this depends on: Phase 3 (#1706) has to
give `Specimen` and `Extraction` pages to hang it off, and the tag → specimen path runs
tag → `node.get_samples()` → `Sample.extraction` → `Specimen`, since `VariantTag` carries no sample. A
node can hold several samples, so the reverse query is a join rather than a lookup. Worth measuring
before deciding whether `VariantTag` wants a denormalised sample or extraction column.

---

## Order of work

1. **A** — #1713's data fixes and namespace filter, then #1714's fields once triage has seen them.
   Independent of everything else, and #1713 corrects a live wrongness in the existing copy consensus
   regardless of whether the rest lands.
2. **C stage 1 + D's analysis launch point + the `TAG_SOMATIC_REPORTABLE` setting** — this is
   sapath#246, and it is the piece that is useful on its own.
3. **C stage 2** — the per-variant screen, initially with the allele-level decision only.
4. **B** — `copy_scope = GENE`, the gene candidate query, and decision 2 on the stage 2 screen.

#444 (multi-variant reporting) starts from stage 1's list: the same triaged set, taken to a report
instead of one classification at a time.

---

## Deferred

- **First-class gene/disease curation object** with periodic review — #1419's schema option, above.
- **Splitting `literature` into gene-level and variant-level** — variantgrid_private#1102. Until then
  `literature` does not travel at gene scope.
- **Per-key provenance for copied values.** Today the whole patch is marked `SubmissionSource.CONSENSUS`
  and the activity page explains it (`classification/templates/classification/activity.html:65`), but
  nothing records *which* record a value came from. Worth doing when the gene-level object arrives,
  since that is when "where did this text come from and when was it reviewed" becomes answerable.
- **Hiding germline-only keys on a somatic form.** `copy_allele_origin` states which keys are
  germline-only, so the form could use it too. That is a larger change to a form many labs rely on, and
  the copy fix does not need it.

## Open questions

| Question | Why it matters |
|---|---|
| Should external-lab candidates be visible by default, or behind a toggle? | Shariant records are the useful ones, but they widen the list |
| Do `amp:level_a`–`d` belong at gene scope? | They are gene + tumour type in practice; if so they need the same human pick, not a scope change |
