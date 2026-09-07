# Tags → Classifications → Multi-variant Report (sample/patient page)

Written by Claude Fable 5 (claude-fable-5), 2026-09-05; revised by Claude Opus 5 (claude-opus-5), 2026-09-07
Status: in progress

Issues: sacgf/variantgrid_sapath#246 (SomaticReportable — easy classifications),
SACGF/variantgrid#444 (Multi-variant Classification + Reporting).

Clickable mockup (agreed design): https://claude.ai/code/artifact/64f391ea-fa7f-496d-a59c-5264d647d28b
— a copy is checked in at `claude/plans/246_classify_report_mockup.html` (open in a browser; it is
self-contained). It shows the agreed layout, labels and interactions; the real page uses the
project's Bootstrap 4 styling, not the mockup's CSS.

## Summary

A new **Classify & Report** tab on the sample page (and patient page, unioned over the patient's
samples) turns outstanding classification tags into classifications, then builds a multi-variant
report from them. Three-stage funnel on one page:

1. **Tag queue** — variants tagged RequiresClassification / SomaticReportable for this case that
   don't yet have a classification. Each row shows the lab's previous classifications for the
   allele (count + conditions).
2. **Classify** — per row, a launcher dialog lists the previous classifications; the scientist
   manually decides whether one applies (conditions are inconsistent, so this is never automated,
   not even highlighted). "Apply to this sample" creates a new classification in the background by
   copying the previous one's consensus, then shows a link to the new record. "New classification —
   full form" creates the record and opens the normal classification editor in a new tab. A
   "Classify all" wizard steps through the queue row by row.
3. **Report** — tick classifications, pick a `ClassificationReportTemplate`, generate a report
   grouped by gene.

Tagging in the analysis is unchanged — still one click, no sample picker. Sample resolution is
automatic where free and otherwise deferred to classification time, where a sample dropdown
already exists.

## Data

### `snpdb.Tag` — one new field

`Tag.allele_origin_bucket` already exists and drives germline (ACMG) vs somatic (AMP) form launch.
The new field marks which tags feed the classify queue (replaces the magic
`settings.TAG_REQUIRES_CLASSIFICATION` string comparison for queue membership; the setting stays
for the retire-on-classify behaviour below).

```python
class Tag(models.Model):  # existing model, snpdb/models/models.py
    id = models.CharField(max_length=50, primary_key=True)
    retired = models.DateTimeField(null=True, blank=True)
    merged_into = models.ForeignKey('self', null=True, blank=True, on_delete=SET_NULL)
    allele_origin_bucket = models.CharField(max_length=1, choices=TAG_ALLELE_ORIGIN_CHOICES,
                                            default=AlleleOriginBucket.UNKNOWN)
    # NEW: tagged variants appear in the sample/patient page classify queue
    requires_classification = models.BooleanField(default=False)
```

### `analysis.VariantTag` — one new field

```python
class VariantTag(GuardianPermissionsAutoInitialSaveMixin, TimeStampedModel):  # existing model
    ...existing fields...
    # NEW: which sample the tagging is about - the study's proband, from the node it was made in.
    # Filled silently at tag time, left null otherwise - never prompted for. Nullable forever.
    sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL)
    # NEW: the classification that satisfied this to-do. The tagging is resolved rather than deleted,
    # so it stays as the record of what was flagged and what it turned into.
    resolved = models.DateTimeField(null=True, blank=True)
    resolved_by = models.ForeignKey(User, null=True, blank=True, on_delete=SET_NULL,
                                    related_name="resolved_variant_tags")
    resolved_classification = models.ForeignKey(Classification, null=True, blank=True, on_delete=SET_NULL)
```

Whose tagging it is comes from the study, not from who carries the variant: a relative who is HET for the
proband's variant does not need their own classification. `AnalysisNode.get_proband_sample()` already answers
this (it is what `AncestorSampleMixin` nodes auto-populate from) and returns None when a node's ancestors
disagree.

### Migrations

1. Schema: add the new fields.
2. Data: set `requires_classification=True` on `settings.TAG_REQUIRES_CLASSIFICATION` and
   `SomaticReportable` where those tags exist; set `SomaticReportable.allele_origin_bucket` to
   SOMATIC if still UNKNOWN.
3. Data: backfill `VariantTag.sample` from the tagged node's `get_proband_sample()`, leaving null
   where the node has none.

No new models. Classifications, report templates and the copy machinery already exist.

## Queue semantics

**Queue membership** (per sample): `VariantTag` rows where `tag.requires_classification`, the tag
is live, and the tagging belongs to this sample — `sample` FK matches, or `sample` is null and the
tag's analysis contains this sample among its `get_samples()`. A null-sample tagging is shown only
to the case's samples that carry the variant, which keeps the row off the pages it cannot be about
without claiming it is theirs — the row's sample stays unset and the dialog's dropdown makes the
scientist choose. Patient page: union over the patient's samples.

**Done** is `VariantTag.is_resolved` — the tagging carries the classification that satisfied it, and
a withdrawn `resolved_classification` puts it back in the queue. The row also shows a classification
the case already has for the allele, so a variant classified elsewhere arrives with the link.

**After classification**, a `requires_classification` tagging is resolved rather than deleted, so it
stays visible everywhere a tagging shows (which is what `SomaticReportable` already did):
- Resolution is automatic when the classification is of the tagging's own sample, or the analysis is
  about one person — `analysis.variant_tag_operations.classification_resolves_tag`.
- Otherwise the queue row renders "✓ Classified" with the link and a **Clear tag** button, so the
  scientist says which person the classification was for. A relative's HET is exactly this case.
- The audit `LogEntry` is now an UPDATE; the analysis audit log reads "Cleared - classified as …".

## Phases

### Phase 1 — models and tag-time sample capture

Migrations above. Then in `set_variant_tag` (analysis/views/views_json.py), fill
`VariantTag.sample` at creation from the node's proband. No UI change; tagging stays one click.

### Phase 2 — the Classify & Report tab

New lazy-loaded tab (the `sample_files_tab` `data-href` pattern) on `view_sample` and
`view_patient`. The view lives in the **analysis** app (home of `VariantTag`; analysis already
imports classification — classification must not import analysis).

Views/URLs:
- `sample_classify_report_tab(sample_id)` / `patient_classify_report_tab(patient_id)` — renders
  queue + classifications panels. Queue rows carry previous-classification summary (count,
  distinct conditions) from the latest published `ClassificationModification` per classification
  of the same allele, restricted to what the user can see (`filter_for_user`), lab's own first.
- Launcher dialog content: the previous classifications listed with curated date, user, condition,
  classification/tier, interpretation-summary snippet. Straight to full-form creation when there
  are none (no dialog).
- "Apply to this sample": one POST to the existing `create_classification` flow —
  `create_classification_object` already accepts `variant_id`, `sample_id`, `lab`, transcript and
  `copy_from_vcm_id`, and applies `ClassificationConsensus(...).consensus_patch`. Add a JSON
  response variant (record id + URL) alongside the existing redirect so the row can update in
  place with the link.
- "New classification — full form": same endpoint without `copy_from_vcm_id`, opened in a new tab
  (record is created at launch, exactly like today's flow, so the queue row links to it
  immediately; evidence is filled in on the classification page).
- "Classify all" wizard: front-end only — steps the launcher dialog through the queue.

### Phase 3 — multi-variant report (#444 beachhead)

- Classifications panel: checkboxes + template `<select>` over `ClassificationReportTemplate`
  (per the #569 discussion: user picks the template rather than per-lab copies).
- New view: `sample_multi_classification_report(sample_id)` taking template name + selected
  `ClassificationModification` ids; groups by gene symbol; renders through the existing report
  template machinery extended to take a list (template context gains `classifications` grouped
  by gene alongside the existing single-record context, so existing single-variant templates
  keep working).
- Structure per #569: case/patient header, summary of results, then gene-grouped variant
  sections.

Linked-classification semantics (compound het linkage, ID linkage discussed on #444) are out of
this plan — selection at report time covers the reporting need without new data.

## Testing

- Proband resolution: a single-sample node names its sample even when that sample doesn't carry the
  variant; a cohort node with one carrier still resolves to nothing.
- Queue query: tag with sample FK; a null-sample tag offered to each carrier with no sample set;
  tag done once resolved, and back in the queue after withdraw.
- Apply-to-sample POST: creates record with consensus patch + sample, returns JSON link.
- Resolve-by-sample keeps the audit LogEntry fields intact, leaves another case's tagging alone, and
  leaves an ambiguous tagging for the Clear tag button.
- URL tests for the two tab views via `URLTestCase`.
