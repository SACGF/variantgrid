# VCF / Cohort page unification and membership editor rework

Written by Claude Fable 5 (claude-fable-5), 2026-09-05

Issue: https://github.com/SACGF/variantgrid_private/issues/933 ("Clean up cohort page" —
which already asked for exactly the chosen editor: autocomplete + a reorderable list with a
remove X per sample). Also folds in variantgrid issues #718 (add whole VCFs to cohorts) and
#1238 (cohort page looks like an error while importing).

Interactive mock-up of the chosen membership editor (variation A "Membership table"):
https://claude.ai/code/artifact/408c3e22-0d3f-407f-a69b-30fd31292744

## Data

**No schema changes.** Every phase works against the existing models; nothing here needs a
migration. The models the plan revolves around, as they stand:

```python
class Cohort(...):
    name = models.TextField()
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)
    version = models.IntegerField(null=False, default=0)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    vcf = models.OneToOneField(VCF, null=True, on_delete=CASCADE)   # NULL for custom multi-VCF cohorts
    parent_cohort = models.ForeignKey("self", null=True, related_name="sub_cohort_set", on_delete=DO_NOTHING)
    sample_count = models.IntegerField(null=True)

class CohortSample(models.Model):
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    cohort_genotype_packed_field_index = models.IntegerField()
    sort_order = models.IntegerField()

    class Meta:
        unique_together = ('cohort', 'cohort_genotype_packed_field_index')
```

The one new data shape is the batch membership payload (replaces per-click
`cohort_op=add/remove` POSTs):

```python
# POST to cohort_sample_edit (or a renamed successor)
{
    "sample_ids": [301, 201, 202, 101],   # full desired membership, in display order
}
```

The response reuses the existing shape from `create_cohort_genotype`
(`status` / `celery_task`) so the client can poll with the standard task machinery.

## Decision record

- **Deferred rebuild stays.** `CohortGenotype` packs per-sample values into fixed-width
  arrays keyed by `cohort_genotype_packed_field_index`, so adding/removing one sample touches
  every row of the partition regardless. An incremental UPDATE path would cost the same I/O as
  a rebuild plus dead-tuple bloat and much trickier code. Membership edits are therefore
  batched in the UI and committed with **one** save that triggers one rebuild (or the free
  sub-cohort conversion). Do not "optimise" this into incremental CohortGenotype updates.
- **One page, two entry URLs.** `view_vcf` stays the canonical page for VCF-backed data (it
  is the only valid address while a VCF is importing / `REQUIRES_USER_INPUT`, or for legacy
  VCFs with no cohort); `view_cohort` keeps redirecting to it when `cohort.vcf` is set and
  renders the same unified template for custom/sub cohorts. Merge happens at the template
  level (tab conditionals + shared includes), with the two views kept as thin functions over
  shared context helpers.
- **Membership editor = variation A** from the mock-up: a single membership table
  (autocomplete add, add-all-from-VCF, drag reorder, strike-through pending removal with
  undo) plus a save bar that shows the pending diff, the single version bump, and — before
  committing — whether Save takes the sub-cohort shortcut or a rebuild.

## Phase 1 — batch membership editing (`views_cohort.py`, `models_cohort.py`)

1. **`Cohort.set_samples(ordered_sample_ids)`**: diff against current membership and apply
   with one version bump.
   - Removals via queryset `.delete()` (bypasses `CohortSample.delete()`'s per-row
     `increment_version()`); additions via `bulk_create` (ditto for `save()`).
   - `sort_order` = position in the list for every row.
   - `cohort_genotype_packed_field_index` for new rows: next free indexes above the current
     max (respecting the unique constraint); `increment_version()` renumbers custom cohorts
     anyway. For a **sub-cohort**, copy the parent's packed index per sample (as
     `create_sub_cohort` does) so the parent CGC keeps being reusable.
   - Exactly one `increment_version()` at the end (which also handles detaching
     `parent_cohort` if a sample outside the parent was added).
2. **`cohort_sample_edit`** accepts the batch payload, calls `set_samples`, then immediately
   calls `create_cohort_genotype_and_launch_task` (which already picks sub-cohort conversion
   vs rebuild) and returns its `status`/`celery_task` — no separate "Save cohort" button for
   the user to forget. Keep permission/archive guards as they are.
3. **Template**: replace the dual-grid tab in `view_cohort.html` with the membership table
   per the mock-up. Client keeps pending state (order array + removed set), computes the
   diff badge and the sub-cohort-vs-rebuild preview (all member samples share one VCF →
   sub-cohort), and submits once. Poll via `wait_for_task`-style JS instead of the inline
   `CC_STATUS_*` state machine. JS goes in a static file under
   `variantgrid/static_files/default_static/js/`, not inline in the template.
4. **Add box**: the existing `SampleChoiceForm` autocomplete, plus an "add all samples from
   VCF…" select (closes the gap in #718). Both only stage client-side; nothing posts until
   Save.
5. **Retire `cohort_sort`** (view, URL, template) — ordering is drag-and-drop in the table
   and part of the same batch payload. Remove `cohort_sample_edit`'s old per-click protocol
   once the new page is in.
6. Import/processing status renders in the save bar as informative progress rather than an
   error-styled message (#1238).

## Phase 2 — unified page from shared includes

1. Extract shared includes used by both pages:
   - `_cohort_samples_table.html` with two modes: VCF mode (editable formset: name, patient,
     extraction, het/hom, detected sex) and cohort mode (read-only rows + membership
     controls from phase 1).
   - `_cohort_wizards.html` (trio/duo/quad wizard block — currently duplicated).
   - The related-data / related-analyses footer block.
2. One outer template, `{% extends base_template %}` with the base (`menu_data_base` vs
   `menu_patients_base`) passed in context. Tab conditionals:
   - always: Details, Samples, Sharing/Permissions (permission target = VCF or Cohort);
   - `{% if vcf %}`: Stats, VCF Info, Relate, Ancestry, Skipped Annotation, downloads,
     reload/archive;
   - `{% if not vcf or cohort.is_sub_cohort %}`: Add/Remove Samples. A sub-cohort's editor
     restricts the add box to the parent VCF's samples (tick/untick — always the instant
     path). A plain VCF cohort states that membership is fixed by the file.
3. **Provenance line** under the page title (extends the existing "This is a sub cohort
   of …" text): one sentence saying where membership comes from, phrased as capability, not
   restriction. VCF cohort: "Samples come from `<vcf>` — membership fixed by the file;
   select samples → Create cohort for a subset." Sub-cohort: "Sub-cohort of `<vcf>` —
   N of M samples; membership edits apply instantly." Custom: "Custom cohort — samples
   drawn from N VCFs; membership edits rebuild genotype data on save." This also explains
   the silent custom→sub-cohort conversion done by `create_cohort_genotype_and_launch_task`.
   Plain text, not an alert box (#1238); the save bar covers the same distinction while
   editing.
4. `view_vcf` and `view_cohort` become thin views sharing context-builder helpers; the
   redirect for VCF-backed cohorts stays.
5. Delete `view_cohort.html`'s bespoke CSS/JS along the way.

## Phase 3 — parity polish and extensions

- Stats tab for custom cohorts: reuse the zygosity graph off the cohort's
  `CohortGenotypeCollection` stats (same plotly helpers as the VCF page).
- **Phenotype search in the add autocomplete** (future extension, user-requested): extend
  the sample autocomplete's search to `patient__phenotype` free text and the matched
  ontology terms behind `PatientTextPhenotype` (annotation.models_phenotype_match), so
  "HP:0001639" or "cardiomyopathy" finds samples whose patients match. Server-side change
  to the autocomplete queryset only; UI unchanged.

## Testing

- `snpdb/tests`: unit tests for `Cohort.set_samples` — single version bump for N-sample
  diff, sort_order matches payload order, packed index copy for sub-cohorts, parent
  detachment when adding an outside sample, empty-list rejected. Fixtures via
  `tests/utils/fake_cohort_data.py:create_fake_cohort`.
- `tests/test_urls.py`: update for the new/removed URLs (`cohort_sort` gone, batch
  `cohort_sample_edit`).
- Existing rebuild/sub-cohort behaviour is covered by `create_cohort_genotype_and_launch_task`
  paths; add a test that a batch edit landing entirely inside one VCF converts to a
  sub-cohort without creating a CGC.
