# #1839 — Analysis tags: one tagging per sample, and the grid says which

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-07; implemented by Claude Opus 5 (claude-opus-5), 2026-09-07
Status: in progress - §2-§5 implemented, under review

[#1839](https://github.com/SACGF/variantgrid/issues/1839)

In an analysis you can see the tags on a variant, but not which sample each tagging is about. That matters when
deciding whether to tag again under the current node: a tagging made in a sibling's branch, or one with no sample at
all, is not this proband's to-do. This plan makes a tagging's sample part of its identity - the same tag can be on a
variant once per sample, plus at most one sample-less tagging - and has the analysis grid draw each pill against the
current node's proband so the difference is visible.

---

## 1. What exists

- `analysis/models/models_variant_tag.py:VariantTag` already carries `sample` - the tagged node's proband at tag time
  (`AnalysisNode.get_proband_sample` of the tagged node), null when the node has no unambiguous proband and
  never prompted for.
- A tagging's identity is the `get_or_create` key in `analysis/views/views_json.py:set_variant_tag`:
  variant, tag, genome build, location, analysis, user. Sample is not in it, so tagging the same variant with the same
  tag from a node with a different proband finds the existing row, moves its `node` and leaves its sample as it was.
  The sample is only ever set on create. "Re-tag here to assign the sample" therefore only works after removing the
  pill first, and then it *moves* the tagging away from whoever it was for.
- Removing a pill deletes every tagging of that tag on the variant in the analysis, all users included.
- The analysis page ships `variantTags = {variant_id: [tag_id]}` and `variantTagsResolved = {variant_id: {tag_id: date}}`
  (`analysis/templatetags/user_tag_color_tags.py:VariantTagsJSNode`, `VariantTagsResolvedJSNode`) and the pill renderer
  `VariantGridFormat.tags` in `variantgrid/static_files/default_static/js/variantgrid_formats.js` draws from those
  alone. The sample is dropped before it reaches the browser. The same dictionaries feed the sample page's variants
  tab (`snpdb/templates/snpdb/data/sample_variants_tab.html`, read only).
- The node grid page knows its node (`var nodeId` in `analysis/templates/analysis/node_data/node_data_grid.html`) but
  not the node's proband.

## 2. Data

### 2.1 `VariantTag` identity

No new field. Sample joins the identity of an analysis tagging, enforced for analysis taggings only - the global
(variant page) taggings on this box already have thousands of same-key duplicates from `analysis/fake_variant_tags.py`
and stay as they are.

```python
class VariantTag(GuardianPermissionsAutoInitialSaveMixin, TimeStampedModel):
    variant = models.ForeignKey(Variant, on_delete=PROTECT)
    tag = models.ForeignKey(Tag, on_delete=CASCADE)
    analysis = models.ForeignKey(Analysis, null=True, on_delete=SET_NULL)
    user = models.ForeignKey(User, on_delete=CASCADE)
    sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL)
    ...  # unchanged fields

    class Meta:
        constraints = [
            # One tagging per (variant, tag, analysis, user, sample); a null sample counts as a value, so
            # there is at most one sample-less tagging as well. Global taggings are outside this.
            models.UniqueConstraint(fields=["variant", "tag", "analysis", "user", "sample"],
                                    nulls_distinct=False,
                                    condition=Q(analysis__isnull=False),
                                    name="varianttag_one_per_sample_in_analysis"),
        ]
```

`nulls_distinct=False` needs Django 5+ and PostgreSQL 15+; this box is Django 6.1 on PostgreSQL 16. The migration adds
the constraint only - it was verified on vg-test2 that no analysis tagging violates it (4 analysis taggings, none with
a sample). Deployments with older data hold repeats a tag merge left behind, and a `ManualOperation` cannot clear them
first - the migrator runs `migrate` before it surfaces manual tasks - so the migration deletes them itself before the
`AddConstraint`, keeping the earliest of each set (the `variant_tags delete-duplicates` rule).

### 2.2 What the browser holds

`variantTags` changes shape from a list of tag ids per variant to a list of taggings, and absorbs the resolved
dictionary - resolution is per tagging, so with one pill per tagging the "done only when every tagging of the tag is
done" rule in `VariantTagsResolvedJSNode` is no longer needed.

```js
// analysis.html / sample_variants_tab.html
variantTags = {
    variant_id: [
        {id: 123, tag: "Consider", sample: 45, resolved: "2026-09-01"},   // resolved: null while to-do
    ],
};
analysisSamples = {45: "Proband", 46: "Mother"};   // id -> name, for tooltips
```

```js
// node_data_grid.html, set per grid render alongside nodeId
var nodeProbandSampleId = 45;   // null when the node has no unambiguous proband
```

### 2.3 Tag endpoint payload

`set_variant_tag` on add always returns the tagging the click landed on, created or not, so the browser can add a pill
exactly when a row was made:

```json
{"variant_tag": {"id": 123, "tag": "Consider", "sample": 45, "sample_name": "Proband", "resolved": null},
 "created": true,
 "node_count_types": [...]}
```

On delete it takes `variant_tag_id` (the pill's own row) instead of variant + tag.

## 3. Rules

**Identity.** Adding tag T to variant V from node N looks for the tagging (V, T, analysis, user, sample = N's proband).
Found: stamp `node` / `node_version` / `node_live_data_sources` as now, and nothing else changes. Not found: create it
with that sample. A node with no proband finds or makes the sample-less tagging. A tagging never changes sample - a
second sample gets a second tagging - so tagging for the proband never takes the tag away from a sibling, and the
sample-less tagging stays as "no one's yet".

**Pills.** One pill per tagging. Its state is decided against the node the grid is showing:

| tagging.sample vs node proband | look | tooltip |
|---|---|---|
| equal (both null included) | as today | `Tagged as T` (+ ` - classified <date>` when resolved) |
| a different sample | small person marker on the pill | `Tagged as T for <sample name>` |
| null, node has a proband | small hollow-person marker | `Tagged as T, no sample - tag here to make one for <proband>` |

The common case looks unchanged. Two pills of the same tag side by side (proband + sibling, or proband + sample-less)
are told apart by the marker. Resolved stays the faded pill with the tick, per tagging.

**Remove.** The X on a pill deletes that one tagging by pk. The analysis' write permission still governs, as now.

**Classify queue.** `analysis/classify_report.py:ClassifyReportCase.variant_tags` shows a case's own-sample taggings
plus sample-less ones the case's samples carry. Where the case already has a tagging of the same variant, tag and
analysis for one of its samples, the sample-less one is left out of that case - the precise one supersedes it there. It
still shows for the other samples of the analysis, which is what "no one's yet" means.

**Resolution.** Unchanged: `analysis/variant_tag_operations.py:classification_resolves_tag` already resolves the
tagging whose sample the classification is of, and a sample-less tagging only in a one-sample analysis.

**Export.** `analysis/tasks/analysis_grid_export_tasks.py` joins tag ids per variant; it dedupes so a tag on two
samples exports once. `analysis/grid_export.py:_summarise_tags_global` is untouched.

**Counts.** Tag node membership (`analysis/models/nodes/filters/tag_node.py:TagNode.tagged_variants_q`) and node counts
are by variant, so a second tagging of the same variant changes nothing there. The toolbar's number of tagged variants
(`setNumVariantTags` in `variantgrid/static_files/default_static/js/analysis.js`) keeps counting keys of `variantTags`.

## 4. Steps

1. **Model** - `Meta.constraints` on `VariantTag` as §2.1, migration `analysis/migrations/`. Run `scripts/vg map`.
2. **Endpoint** - `analysis/views/views_json.py:set_variant_tag`: the sample is worked out from the node *before* the
   lookup (`get_proband_sample`) and is part of the key; return the §2.3 payload. Delete by `variant_tag_id` from the
   analysis path, checking the tagging belongs to the analysis. The variant page path (no analysis) is unchanged.
3. **Page state** - `analysis/templatetags/user_tag_color_tags.py`: `VariantTagsJSNode` emits §2.2 in one query
   (`values_list` of id, variant, tag, sample, resolved, resolved_classification__withdrawn - resolved is null in the
   payload where `VariantTag.unresolved_q` would keep it a to-do). Delete `VariantTagsResolvedJSNode` and
   `render_variant_tags_resolved_dict`; both templates drop the second line. Add `render_analysis_samples_dict` for
   `analysisSamples` from `analysis/models/models_analysis.py:Analysis.get_samples`.
4. **Node grid** - `analysis/views/views_node.py:node_data_grid` puts `node_proband_sample_id` in the context;
   `node_data_grid.html` sets `nodeProbandSampleId` next to `nodeId`. The sample variants tab sets it to the page's
   sample, so its read-only pills are drawn against that sample.
5. **Renderer** - `VariantGridFormat.tags`: one pill per entry, state per §3, pk on the pill (`variant_tag_id`
   attribute) for delete. `getVariantTagHtml` in `variantgrid/static_files/default_static/js/grid.js` grows the marker
   and pk. `sortVariantTags` sorts entries by their tag. `tagClickHandler` / `removeVariantTag` send the pk.
6. **Add / remove state** - `setVariantTag` in `analysis.js`: on add, push the returned tagging only when `created`
   and record its sample name in `analysisSamples`; on delete, drop the entry by pk. The done-marker handling for
   re-tag goes away with `variantTagsResolved`.
7. **CSS** - the marker in `variantgrid/static_files/default_static/css/global.scss` next to `.grid-tag-resolved`,
   mirrored by hand into `global.css`.
8. **Queue** - `ClassifyReportCase.variant_tags` dedupe per §3.
9. **Export** - dedupe per §3.
10. Docs: the tagging identity and the pill rule are two lines in `analysis/CLAUDE.md` next to the existing tag notes.

## 5. Tests

`analysis/tests/test_variant_tags.py` already covers resolution and the resolved dictionary; the latter's test moves
to the new dictionary. Worth keeping:

- Tagging from a node whose proband differs from an existing tagging's sample makes a second tagging; from a node
  with the same proband finds it and stamps the node; from a node with no proband finds or makes the sample-less one.
  (`analysis/tests/test_variant_tags.py`, through `set_variant_tag`.)
- The page dictionary emits one entry per tagging with its sample, and `resolved` null for a withdrawn resolving
  classification.
- The queue leaves the sample-less tagging out of a case that has its own for the same variant and tag, and keeps
  it for a case that does not.

The pill states are JS; check them with `vg page` on an analysis with a trio and by eye, not with a unit test.

## 6. Out of scope

- The variant page tag list and the tag node's Classifications tab list a row per tagging already; a sample column
  there is a separate small change.
- A tagging changing its sample. The model is remove and re-tag.
