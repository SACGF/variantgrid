# Reusing prior curation for somatic reporting

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-08-11; revised by Claude Opus 5 (claude-opus-5),
2026-09-07 and Claude Opus 5.5 (claude-opus-5-5), 2026-10-12
Status: in progress (A, B, C stage 1's list, the stage 2 wizard walk and D's specimen / extraction launch
points landed; stage 1's drop and external-lab marking, stage 2's condition step and bucket seeding, and
D's analysis-scoped button remain)

Design for [#1419](https://github.com/SACGF/variantgrid/issues/1419) (separate out gene / disease
curation), [sapath#246](https://github.com/SACGF/variantgrid_sapath/issues/246) (SomaticReportable →
easy classifications), and the entry point [#444](https://github.com/SACGF/variantgrid/issues/444)
will consume.

The driving observation: in somatic work the same variants recur constantly, and the lab has usually
curated them before. Copy consensus already reuses prior work per allele. What is missing is a way to
*see* what has been curated before you start, a way to get through a report's worth of variants without
re-navigating for each one, and a way to reuse the gene-level content that is identical across every
variant in that gene.

---

## What already exists

Most of the machinery is here; the delta is smaller than the issue titles suggest.

- **Allele-level copy consensus.** `classification/models/classification.py:ClassificationConsensus` returns
  the latest record for an allele per bucket — published, non-withdrawn, `exclude_external_labs=True` —
  defaulting the radio to whichever matches the user's `allele_origin_focus`. `consensus_patch` copies
  `value` and `note` for every ekey whose `copy_scope` is in scope, and forces `allele_origin` from the
  source's bucket rather than copying it.

- **Somatic vs germline is data, not code.** The `allele_origin` option `somatic` carries
  `namespaces: ["somatic"]` and `assertion_method` option `amp` carries `["amp"]`;
  `_evidence_key_overrides_from_evidence_fields` (`classification/models/classification.py:1965`) turns
  those into the per-record namespace set that switches the `amp:`/`horak:`/`somatic:` keys on. So
  "launch as AMP/somatic" means seeding those two keys at creation.

- **Tag → classification.** `analysis/views/views.py:CreateClassificationForVariantTagView` already
  narrows the sample dropdown to the tag's node samples and posts to `create_classification_for_analysis`,
  which links the new record to the analysis via `AnalysisClassification`. Which tags ask for a
  classification is `snpdb/models/models.py:Tag.requires_classification`, set for `RequiresClassification`
  and `SomaticReportable` by `snpdb/migrations/0251_one_off_tags_requiring_classification.py`; every "New
  Classification" button reads that flag rather than naming a tag.

- **Gene lookup.** `ResolvedVariantInfo.gene_symbol` is a real FK
  (`classification/models/classification_variant_info_models.py:143`) and
  `classification_gene_symbol_filter` exists, so "prior classifications in this gene" is one query with
  no new denormalisation.

On the specimen side, `analysis/models/models_variant_tag.py:VariantTag` carries a `sample` (the
study's proband at tag time, sapath#246). Where it is unset — a tagging that predates the field, or a node
with no proband — tag → specimen goes tag → `node.get_samples()` → `Sample.extraction` → `Specimen`.

---

## Part A — which evidence keys copy

`classification/models/evidence_key.py:EvidenceKey` carries `copy_scope` (`NONE` / `ALLELE` / `GENE`) and
`copy_allele_origin` (`ANY` / `GERMLINE` / `SOMATIC`), both in the admin fieldset beside `max_share_level`
(`classification/admin/classification_admin.py`); `classification/management/commands/legacy_somatic.py`
maps its CSV's `copy_consensus` column onto `copy_scope`. The values below are set by data migrations from
two issues:

- **[#1713](https://github.com/SACGF/variantgrid/issues/1713)** — A1 and A2,
  `0170_copy_consensus_off_patient_and_report_keys` and `0171_copy_consensus_off_annotation_backed_keys`.
- **[#1714](https://github.com/SACGF/variantgrid/issues/1714)** — the two fields, A3's twenty
  germline-only keys and Part B's gene scope, `0176_evidence_key_copy_scope` and
  `0177_evidence_key_copy_scope_and_allele_origin_values`.
  `consensus_patch` filters on the source record's allele origin bucket, which is the bucket the new
  record inherits; Part C's wizard seeds the bucket itself, so it wants the target bucket passed in.

A3's namespace filter is not done. Filtering the patch by the target record's namespaces drops all 28
`acmg:*` criteria, because the create form never sets `assertion_method` and the "assume ACMG unless
Horak" rule exists only in `variantgrid/static_files/default_static/js/vc_keys.js`.

### A1 — Test-level facts about *this* tumour: `copy_scope = NONE`

| Key | Category | Why it is not copied |
|---|---|---|
| `somatic:tmb_value`, `somatic:tmb_status` | HT | TMB of the specimen being reported |
| `somatic:msi_value`, `somatic:msi_status` | HT | MSI of the specimen being reported |
| `somatic:hrd_status` | HT | HRD of the specimen being reported |
| `testing_context` | HT | Per-test |

These describe the patient's tumour, not the allele: copying them from a prior classification would import
another patient's measurements. `somatic:tumor_cellularity` is `NONE` for the same reason.

`SpecimenMeasure` (#1559) is the source for the five `somatic:` ones — they arrive through
`evidence_from_sample_and_patient.py`, not through a copy.

### A2 — Report and reviewer free text: `copy_scope = NONE`

`somatic:summary_interpretation` is described in its own migration
(`classification/migrations/0164_somatic_hrd_msi_tmb_ekeys.py`) as a free-text summary of all
classifications in a report, `max_share_level='lab'`, flagged there as very likely to contain
patient-identifiable information. Copying it would pull another patient's report narrative into this
record.

`review_comment` (SO) has the same shape — a reviewer's comment about a different record.

### A3 — Cross-bucket leakage

`consensus_patch` uses `EvidenceKeyMap.cached()`, the unconfigured global map, so it copies every copyable
key regardless of the target record's namespaces. `exclude_namespace` (`evidence_key.py:567`) only
suppresses mandatory-field validation (`classification.py:1241`) — it does not stop a patch writing the
key. So copying a germline record into a somatic one writes 28 `acmg:*` criteria into evidence that the
form will never show, and the reverse writes `amp:`/`horak:` into a germline record.

**Fix: filter the consensus patch by the target classification's `evidence_key_overrides.namespaces`.**
That handles the 28 `acmg`, 17 `horak` and 4 `amp` copyable keys with no per-key data entry, and it stays
scoped to the UI-initiated copy — API imports keep writing out-of-namespace keys deliberately, which is
how records shared from labs with a different config survive.

Germline concepts that have no namespace to filter on carry `copy_allele_origin = GERMLINE`:

```
condition_incidence  mode_of_inheritance  gene_penetrance  variant_penetrance  proband_count
segregation  segregation_meioses  segregation_affectedcarriers  segregation_unaffectedcarriers
segregation_affectednoncarriers  segregation_bayes  segregation_lod  s_other  s_summary
denovo_points  d_other  d_summary  a_other  a_summary  match_maker_exchange
```

Twenty rows of data entry rather than a rule, because there is nothing in the schema that distinguishes
them.

### A4 — Redundant but harmless, left alone

The gnomAD count keys and about a dozen predictor scores (`alphamissense`, `bayesdel`, `vest`,
`varity_r`, `mpc`, `primateai`, `metarnn`, `clinpred`, `mutpred2`, `aloft`) are copyable while their
siblings (`revel`, `cadd`, `sift`, `polyphen`) are not. Inconsistent, but autopopulate runs first and
wins: `create_classification_object` applies the consensus patch with `leave_existing_values=True`
(`classification/views/views.py:428-440`) and `AutopopulateView` guards on `used_keys` (`:349-356`), so
a copied value only ever fills a gap annotation left. Not worth churning.

### The `condition` key

`condition` is copyable. For somatic it holds the tumour type, which is a per-patient fact and the
single axis on which reuse is most likely to go clinically wrong. It stays copyable, but the wizard
surfaces it as an explicit human decision rather than a silent pre-fill (Part C, decision 3).

---

## Part B — gene-level reuse, with the human picking

The gene-level content of a somatic classification — what the gene does, its role in cancer, the
gene-level literature — is identical across every variant in that gene, and is exactly what the
allele-keyed copy cannot reach.

Part A provides the vocabulary (`copy_scope = GENE` on the `H` keys below). This part consumes it: the same
candidate list in three places, each with a different amount of ceremony.

**Landed.** The code is `classification/models/classification.py:ClassificationConsensus` (the candidate
functions and `apply_to`), `classification/views/views.py:CreateClassificationForVariantView` plus
`classification/templates/classification/create_classification_for_variant.html`,
`analysis/classify_report.py` plus `analysis/templates/analysis/classify_report_tag_dialog.html`, and
`classification/views/views_gene_consensus.py` plus the Gene Content card in
`classification/templates/classification/classification.html`. The sections below are the rules those
implement — Part C follows the same ones.

### The bucket is an input, never an output

A germline record is never copied into a somatic one, in either scope. `consensus_patch` takes
`allele_origin` from the *source*, so an unfiltered candidate list would let a somatic-bucket tag copy a
germline record and come out germline.

The rule everywhere: decide the target bucket first, then list only candidates in that bucket.

- **Create-from-variant page** — the bucket comes from the user's `allele_origin_focus`
  (`ClassificationConsensus.default_allele_origin_bucket`), shown as an explicit germline / somatic choice
  at the top of the copy section so it can be flipped. Both buckets' candidates are rendered and the
  inactive set is disabled, so only the chosen bucket's radios can be posted.
- **Classify & Report** — the bucket comes from the tag
  (`analysis/classify_report.py:tag_allele_origin_bucket`; `Tag.allele_origin_bucket` "Both" means no
  filter).
- **In-form helper** — the bucket is the record's own `allele_origin_bucket`.

`consensus_patch` seeds `allele_origin` from the source, which is always the same bucket.
This also takes most of the weight off A3: same-bucket copies mostly share namespaces, so the remaining
cross-namespace leak is `amp` vs `horak` within somatic, and `copy_allele_origin = GERMLINE` covers
germline-only keys that have drifted into somatic records.

Part C's wizard seeds the bucket itself, so it wants the target bucket passed *in* rather than derived —
`ClassificationConsensus.allele_origin_bucket` reads it off the source record today, and that is the one
place to change when C lands.

### The candidates are shown, never auto-picked

AMP tiering and therapy content are gene **and tumour type**, not gene. Copying gene-level content from
a colorectal case onto a melanoma case is a clinical error rather than staleness, and the phenotype data
available to match on is not good enough to automate the judgement. So the gene-level source is always a
human choice, and "none" is the default.

The query: `ClassificationModification.latest_for_user(user, published=True, exclude_external_labs=True)`
filtered through `classification_gene_symbol_filter(gene_symbol)` and the target bucket, excluding any
record already offered as an allele-level candidate. It lives in
`classification/models/classification.py:ClassificationConsensus.gene_consensus_groups`, which all three
surfaces call.

### Deduplicated — one row per distinct gene content

Most records in a gene carry identical gene content, because they were copied from each other. Listing
them all is noise. Group the candidates by the values of their `GENE`-scope keys (`value` and `note`,
after the `copy_allele_origin` filter, so two records that differ only in a variant-level field are the
same group). Each group is one row:

- **Representative** — the most recently curated record in the group. The row shows its curated date,
  clinical significance / tier, condition and lab, which is enough to judge tumour-type relevance.
- **"and N other records"** — a disclosure on the row. Expanded, it lists the rest of the group as one
  overview line each (tier, condition, date, lab), so the curator can see the spread of tumour types this
  gene content has been used for. That spread is often the deciding information.

Groups are ordered by the representative's curated date, newest first, capped at ten groups. Whichever
record is picked, only its `GENE`-scope keys travel (`copy_scopes=COPY_SCOPES_GENE`). A record whose
gene-scope keys are all empty is not a candidate — a row with nothing to copy is noise.

### Where the pick is made

**1. Create-from-variant page** (`classification/views/views.py:CreateClassificationForVariantView`, the
primary place). A second radio group, "Gene information from", under "Copy values from", with the
deduplicated rows above and "none" selected. It posts `copy_gene_from_vcm_id`. The autofill preview labels
each key's source as annotation, allele copy (`copy from latest`) or gene copy (`copy from gene`);
`used_keys` ordering in `AutopopulateView` gives autopopulate over allele over gene.

Picking an allele-level candidate hides the gene group: the allele record's gene content travels with it
(`GENE` is copyable at allele scope already), and there is no mixing and matching between two source
records. The gene group appears when the allele choice is "none", which is the case gene copy exists
for — a variant the lab has never seen in a gene it curates often.

**2. Classify & Report dialog** (`analysis/templates/analysis/classify_report_tag_dialog.html`). "Apply to
this sample" on an allele candidate keeps working exactly as it does and brings that record's gene content
with it. When the tag has no copyable allele-level candidate, the dialog shows the gene rows instead, with
the same button. Both surfaces call `ClassificationConsensus.gene_consensus_groups`, so they cannot drift.

**3. In-form helper box** (`classification/views/views_gene_consensus.py`, rendered as the Gene Content
card in `classification/templates/classification/classification.html`, right column beside Criteria
Summary). For the records the create page never sees: created from the Classify & Report tab with no
candidate, created by API or import, or created before the gene was curated. The card reads "Gene content
available from N classifications in *GENE*" and opens a dialog with the same deduplicated rows, each row
also showing its `GENE`-scope values beside this record's current ones. Apply patches empty fields only,
as `SubmissionSource.CONSENSUS`, then reloads the form. Both card and dialog are AJAX-loaded, so the
candidate query stays off the form's own page load.

The box is present only while it is useful. It is hidden once gene content has been applied to this
record — by the box itself, or by the create page's copy — and the fact is derivable without new data:
a `ClassificationModification` with `source = CONSENSUS` whose `delta` touches a `GENE`-scope key. While
that modification exists the box collapses to one line only when a candidate's curated date is later
than it — "newer gene content from *lab*, *date*" — so a record can learn that the gene has been
re-curated since, and otherwise takes no space at all.

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

Show them, do not copy from them. `all_consensus_candidates` passes `exclude_external_labs=True`, which
hides what other labs — including Shariant — have curated on the same allele, and that is worth seeing
when deciding how to curate. But an external record was curated under another lab's config, assertion
method and namespaces, so copying its evidence into a local record imports assumptions that were never
reviewed here. So the create page has a read-only "Other labs" row
(`ClassificationConsensus.external_lab_candidates`), and the Classify & Report dialog — which lists every
visible lab through `_previous_by_tag` — shows an external record without the "Apply to this sample"
button (`PreviousClassification.can_copy`).

---

## Part C — the wizard

A report's worth of somatic variants is a batch, and classifying it one tag at a time is N independent
navigations. But each variant still needs real decisions, so batching the *decisions* would trade
navigation cost for a screen nobody can reason about. The wizard batches the triage and then serialises
the decisions, one variant at a time.

### Stage 1 — triage the list

**Partly landed.** The list is the Classify & Report tab (`analysis/classify_report.py`,
`analysis/views/views_classify_report.py`), on the sample, patient, specimen and extraction pages: one row
per queue tagging in scope, with the variant, its gene, the tag, the sample and what already exists for
that allele. Two things in this section are still open:

- **The drop gesture** below — a row has no control to remove a tagging the lab will not report.
- **The row's overview** shows only the count of previous classifications and their conditions. It still
  needs the latest one's clinical significance / tier, lab and date, with external labs marked as such
  (the dialog already marks them, "another lab – for context").

A row leaves the queue by being *resolved* against a classification
(`analysis/variant_tag_operations.py:resolve_variant_tag`) rather than deleted, so the tagging stays as
the record of what was flagged. That is the "classified" exit; the "not reporting this" exit is the drop
below, which does not exist yet.

One row per unresolved tagging (`analysis/models/models_variant_tag.py:VariantTag.unresolved_q`) in scope
whose tag is a somatic classify-queue tag — `snpdb/models/models.py:Tag.classify_queue_qs_for_bucket`
with `AlleleOriginBucket.SOMATIC`, so a tag marked "Both" counts too. Each row shows the variant and its
gene, and a brief overview of what already exists for that allele: how many prior classifications, the
latest one's clinical significance / tier, its condition, lab and date — external labs included, marked
as such.

The scientist drops the rows they are not going to report. Dropping deletes the tagging, which is the
same gesture as untagging in the analysis — a plain delete, with the existing post-delete signal
(`analysis/signals/signal_handlers.py:29`) updating node counts and nothing else recorded. Then:
*Classify N variants*.

The drop is a hard delete because a soft-delete flag on `VariantTag` would have to be honoured by 32 query
sites across 15 modules (grids, node counts, `get_for_build`, liftover tasks, the analysis JS tag dict), and
`get_or_create` in `set_variant_tag` (`analysis/views/views_json.py:227`) would need to undelete rather
than duplicate.

Classifying resolves the tagging exactly as the germline queue does
(`analysis/variant_tag_operations.py:resolve_variant_tag`): the row is stamped rather than deleted, so
the tag survives the classification that came from it — which is what sapath#246 asks for, the somatic
lab keeping its tags on the variants it has reported.

### Stage 2 — one variant at a time

**Partly landed.** "Classify all (wizard)" on the tab walks the unclassified rows in order, one
`analysis/templates/analysis/classify_report_tag_dialog.html` dialog each, with a position ("tag 4 of 11")
and Skip; a skipped row is not offered again in the same walk
(`variantgrid/static_files/default_static/js/analysis/classify_report_tab.js`). The dialog carries
decisions 1 and 2 below as "Apply to this sample" buttons — allele candidates, or gene groups when the
allele has none — and posts to `analysis/views/views_classify_report.py:create_classification_for_case`.
Sample and lab are shown as plain fields. Still open:

- **Decision 3, the condition** — the dialog has no condition field.
- **Seeding the bucket** — `create_classification_for_case` seeds neither `allele_origin` nor
  `assertion_method`, so a "none" pick creates a record with no bucket. This is where the target bucket
  gets passed into `ClassificationConsensus` (Part B, "The bucket is an input").

The screen as designed:

1. **Which prior classification for this allele to copy from** — radio over the candidates, each showing
   condition, clinical significance / tier, lab, date. Plus "none".
2. **Which prior classification for this gene to copy from** — same shape, only shown when the gene has
   candidates, defaulting to none (Part B).
3. **The condition for this record** — pre-filled from whichever record was picked in 1 or 2, always
   visible and always editable.

Lab, sample and transcript resolve silently: lab from user settings, sample from the tag's node
(existing `_get_sample_form` behaviour, `analysis/views/views.py:1223`), transcript from the canonical
selection. The dialog shows the transcript and lets sample and lab be changed for the cases where the
default is wrong. Create, and the wizard advances to the next variant.

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

**The tag vocabulary** is `Tag.requires_classification` and `Tag.allele_origin_bucket` on each tag — a
deployment flags its own tags on the tag settings page, and a lab can have several. A deployment turns
this feature on by flagging a tag as a somatic queue tag there; SA Path's are, from `snpdb/migrations/0251_one_off_tags_requiring_classification.py`. Every
launch point counts and scopes with the same `Tag` classmethods as stage 1's list.

**Landed launch points.** The Classify & Report tab on the sample, patient, specimen and extraction pages
(`analysis/views/views_classify_report.py`, case types `sample` / `patient` / `specimen` / `extraction`;
`patients/templates/patients/view_specimen.html`, `patients/templates/patients/view_extraction.html`). A
specimen or extraction case resolves to its samples through `patients/sample_grouping.py:get_sample_group`.
From the analysis, each sample pill in the Input Samples panel carries a classification count linking to
that sample's tab (`analysis/templates/analysis/analysis_input_samples.html`).

Still to do:

**From the analysis**, beside the existing tags button — "Classify somatic (N)" — scoped to that
analysis's taggings rather than the whole sample's.

---

## Remaining work

1. Stage 2's condition field and `allele_origin` / `assertion_method` seeding, with the target bucket
   passed into `ClassificationConsensus`.
2. Stage 1's drop, and the row's tier / lab / date with external labs marked.
3. D's analysis-scoped "Classify somatic (N)" button.

#444 (multi-variant reporting) builds from the same tab: its "Build multiple variant report" card sits
under stage 1's list.

---

## Deferred

- **First-class gene/disease curation object** with periodic review — #1419's schema option, and the
  only thing that can express "reviewed on this date, due for review". It needs the disease axis settled
  first (condition matching, ontology hierarchy, whether a gene/disease record is per-lab). The copy route
  does not foreclose it: `copy_scope = GENE` is exactly the set of keys such an object would own.
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
| Do `amp:level_a`–`d` belong at gene scope? | They are gene + tumour type in practice; if so they need the same human pick, not a scope change |
