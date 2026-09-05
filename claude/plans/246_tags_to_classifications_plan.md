# Tags → Classifications → Multi-variant Report (sample/patient page)

Written by Claude Fable 5 (claude-fable-5), 2026-09-05

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
    # NEW: which sample the tagging is about. Filled silently at tag time when unambiguous,
    # left null otherwise - never prompted for. Nullable forever (old/ambiguous tags).
    sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL)
```

### Migrations

1. Schema: add both fields.
2. Data: set `requires_classification=True` on `settings.TAG_REQUIRES_CLASSIFICATION` and
   `SomaticReportable` where those tags exist; set `SomaticReportable.allele_origin_bucket` to
   SOMATIC if still UNKNOWN.
3. Data: backfill `VariantTag.sample` — rule (a) the tag's node has exactly one sample →
   that sample; rule (b) otherwise, exactly one of the analysis's samples has the variant
   (zygosity not ref/missing) → that sample; else leave null.

No new models. Classifications, report templates and the copy machinery already exist.

## Queue semantics

**Queue membership** (per sample): `VariantTag` rows where `tag.requires_classification`, the tag
is live, and the tag resolves to this sample — `sample` FK matches, or `sample` is null and the
tag's analysis contains this sample among its `get_samples()` and the sample carries the variant.
Patient page: union over the patient's samples.

**Done** is inferred, never stored on the tag: a classification exists for the tag's allele with
`Classification.sample` among the case's samples. A null-sample (ambiguous) tag counts as done
once *any* of its analysis's samples has a matching classification. Withdrawn classifications
don't count, so the tag reappears if one is withdrawn.

**After classification**, existing behaviour per tag is kept:
- RequiresClassification is a to-do — retired (deleted with audit log) via
  `analysis.variant_tag_operations.retire_requires_classification_tags`. That function filters by
  analysis; extend it (or add a sibling) to retire by variant + case samples when triggered from
  the sample page.
- SomaticReportable stays (somatic labs use it for tracking) — its queue row renders as
  "✓ Classified" with a link to the classification.

## Phases

### Phase 1 — models and tag-time sample capture

Migrations above. Then in `set_variant_tag` (analysis/views/views_json.py), fill
`VariantTag.sample` at creation using backfill rules (a)/(b). No UI change; tagging stays
one click.

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

- Backfill rules (a)/(b) and null outcome on a cohort with two carriers.
- Queue query: tag with sample FK, tag resolved by inference, tag excluded once a matching
  classification exists, tag reappearing after withdraw.
- Apply-to-sample POST: creates record with consensus patch + sample, returns JSON link.
- Retire-by-sample extension of `retire_requires_classification_tags` keeps the audit LogEntry
  fields intact.
- URL tests for the two tab views via `URLTestCase`.
