# #21 — "Tagged" node count

[#21](https://github.com/SACGF/variantgrid/issues/21) asks for tag counts on analysis nodes. This plan
delivers a single **Tagged** count - the number of variants in a node that carry any tag *in this
analysis* - as a new built-in node count that users switch on in the same place as ClinVar / OMIM /
Classified. Per-tag counts can be layered on later through the same settings page once the analysis-scoped
count and its bookkeeping exist.

## How node counts work today

- `BuiltInFilters` (`snpdb/models/models_enums.py:79`) is the fixed vocabulary: single-letter codes,
  `CHOICES` (with `TOTAL`) drive the node count settings, `FILTER_CHOICES` (without `TOTAL`) drive
  `BuiltInFilterNode`, `COLORS` drives the legend CSS (`user_tag_color_tags.render_node_count_colors_css`).
- A node's counts are computed in one aggregate inside the Celery load - `get_node_counts_and_labels_dict`
  (`analysis/models/nodes/node_counts.py:107`) builds one `Count(filter=q)` per configured label from
  `get_extra_filters_q` and writes `NodeCount(node_version, label, count)` rows
  (`AnalysisNode.node_counts`, `analysis/models/nodes/analysis_node.py:1046`).
- Labels come from `Analysis.get_node_count_types()` (`analysis/models/models_analysis.py:245`), falling
  back to `BuiltInFilters.DEFAULT_NODE_COUNT_FILTERS`. Enabling a new label in settings reloads every
  node (`changeAnalysisSettings`, `analysis_nodes.js:1015`).
- Clicking a count opens the node grid with `extra_filters=<code>`; `VariantGrid._get_q`
  (`analysis/grids.py:169`) applies the same `get_extra_filters_q` Q.
- The client paints counts from `nodes_status` (`analysis/views/views_json.py:372`) via
  `messagePoller.observe_node(node_id, "count", ...)` inside `updateDirtyNode` (`analysis_nodes.js:593`).
  A node has two client-visible versions: `version` (queryset changed - node reloads) and
  `appearance_version` (draw it again - `retrieveAndUpdateNodeAppearances`, `analysis_nodes.js:1072`).
- Tagging: `set_variant_tag` (`views_json.py:213`) creates/deletes `VariantTag` rows. The `post_save` /
  `post_delete` handlers (`analysis/signals/signal_handlers.py:19`) set visible `TagNode`s dirty
  synchronously and hand the rest to Celery. `VariantTag.variant` is always in the analysis build, and
  `TagNode._get_node_q` (`analysis/models/nodes/filters/tag_node.py:59`) matches analysis tags on
  `variant_id` for that reason.

## Design

Two halves:

1. **The count itself** is an ordinary built-in filter: `Q(pk__in=VariantTag.objects.filter(analysis=...)
   .values("variant_id"))`. It rides the existing aggregate, so every node reload computes it from
   scratch for free. Being a subquery keeps it a subset of the node's rows, which the
   `label counts > total` check at `analysis_node.py:1093` relies on.
2. **Keeping it current between reloads.** A tag add/remove changes the Tagged count of every node that
   holds the variant, and nodes hold their counts in `NodeCount` rows keyed by an immutable
   `node_version`. So tagging adjusts those rows in place (`count = F("count") ± 1`) and bumps
   `appearance_version` on the nodes that changed. The node's `version` is untouched, so nothing
   reloads; the next real reload recomputes the label from scratch, which also heals any drift.

   Membership test per node: nodes at or under `ANALYSIS_NODE_STORE_ID_SIZE_MAX` (1000) already store
   their exact pk set on the TOTAL row (`AnalysisNode.get_cached_node_pks`, `analysis_node.py:642`) -
   a Python `in`. Larger nodes get `node.get_queryset().filter(pk=variant_id).exists()`, a pk-anchored
   probe. Analyses are 10-30 nodes, so a tag click costs a few row updates.

   "Tagged" counts *variants*, and a variant can carry several tags, so the row moves only when the
   variant's tag count in this analysis crosses zero: on add, when this is the variant's first tag in
   the analysis; on delete, when it was the last.

## 1. `BuiltInFilters` and migrations

`snpdb/models/models_enums.py:79`:

```python
TAGGED = 'A'  # Any tag in this analysis

# FILTER_CHOICES stays as-is (BuiltInFilterNode - TagNode already covers tags as a filter)
CHOICES = [(TOTAL, 'Total')] + FILTER_CHOICES + [(TAGGED, 'Tagged')]
COLORS = [..., (TAGGED, "#f5a623")]
```

Any unused colour that reads as distinct from the existing black/red/green/lavender/purple/pink/mauve/blue
works; orange is a suggestion.

`CHOICES` is the `choices=` on `NodeCountSettings.built_in_filter`
(`snpdb/models/models_user_settings.py:703`) and `AnalysisNodeCountConfigRecord` via
`AbstractNodeCountSettings`, so `makemigrations` produces two choices-only migrations
(`snpdb/0224_...`, `analysis/0125_...`). Node view label (`analysis/views/nodes/node_view.py:39`),
settings tabs (`get_node_counts_mine_and_available`) and the legend all pick the new entry up from
`CHOICES` / `COLORS` with no further change.

`DEFAULT_NODE_COUNT_FILTERS` stays `[TOTAL, IMPACT_HIGH_OR_MODERATE, CLINVAR]`: Tagged is opt-in per
analysis (settings cog → Node counts) or as a user default (user settings → node counts).

## 2. `get_extra_filters_q` takes the analysis

`analysis/models/nodes/node_counts.py:26` currently takes `(user, annotation_version, extra_filters)`.
The Tagged Q needs the analysis, and both existing arguments are analysis attributes, so change the
signature to `get_extra_filters_q(analysis, extra_filters)` and derive `user` / `annotation_version`
inside. Add the branch:

```python
elif extra_filters == BuiltInFilters.TAGGED:
    q = Q(pk__in=Variant.objects.filter(varianttag__analysis=analysis).values("pk"))
```

Written against `Variant` (from `snpdb.models`) via the reverse `varianttag` relation because
`models_variant_tag` imports `analysis_node`, which imports `node_counts` - importing `VariantTag` here
would cycle. The subquery is analysis-scoped so it's small, and `pk__in` needs no `distinct()`.

Update the three callers: `get_node_counts_and_labels_dict` (`node_counts.py:113`, has `node.analysis`),
`BuiltInFilterNode._get_node_q` (`analysis/models/nodes/filters/built_in_filter_node.py:30`),
`VariantGrid._get_q` (`analysis/grids.py:172`).

## 3. Cached-count shortcut

`AnalysisNode._get_cached_label_count` (`analysis_node.py:996`) reuses the `cloned_from` node version's
counts. `Analysis.clone()` (`models_analysis.py:285`) copies nodes with `cloned_from` set but leaves
`VariantTag`s behind, so the source analysis's Tagged count is meaningless for the copy. Skip the
`cloned_from` branch when `label == BuiltInFilters.TAGGED`. The parent-derived shortcuts (all parents
zero → zero; unmodified single parent → parent's count) stay valid, since parent and child share the
analysis.

## 4. Adjusting counts on tag add / remove

New module `analysis/variant_tag_node_counts.py` (sits beside `variant_tag_operations.py`; it needs
`Analysis`, `AnalysisNode`, `NodeCount`, `VariantTag` and `BuiltInFilters`):

```python
def update_tagged_node_counts(analysis: Analysis, variant_id: int, delta: int):
    """ A variant gained its first / lost its last tag in this analysis - move the Tagged count on every
        loaded node that holds it, without a reload """
```

- Early exit when `BuiltInFilters.TAGGED` is absent from `analysis.get_node_count_types()` - analyses
  without the count pay nothing.
- Iterate `analysis.analysisnode_set.select_subclasses()`. For each node fetch its current
  `NodeCount(node_version=node.node_version, label=TAGGED)`; a node with no such row (loading, errored,
  freshly dirtied `TagNode`, count added after the last load) is skipped - its next load computes the
  label from scratch.
- Membership: `AnalysisNode.get_cached_node_pks(node)` → `variant_id in pks`; otherwise
  `node.get_queryset().filter(pk=variant_id).exists()`.
- For holders: `NodeCount.objects.filter(pk=...).update(count=F("count") + delta)` and
  `AnalysisNode.objects.filter(pk=node.pk).update(appearance_version=F("appearance_version") + 1)` -
  queryset updates so no `save()` side effects (dirty flags, child propagation) fire.

Call sites in `analysis/signals/signal_handlers.py`, after `analysis_tag_nodes_set_dirty(..., visible=True)`
so dirtied `TagNode`s are already on a new version and fall through the "no row" skip:

- `variant_tag_create`: when `instance.analysis` and no other `VariantTag` for
  `(analysis, variant_id)` exists → `delta=+1`.
- `variant_tag_delete`: when `instance.analysis` and no `VariantTag` for `(analysis, variant_id)`
  remains (handler is `post_delete`, so the row is already gone) → `delta=-1`.

Both run synchronously in the request, like the visible-`TagNode` dirtying beside them and for the same
reason: the client calls `checkAndMarkDirtyNodes` in the tag success callback (`analysis.js:748`) and
needs the bumped `appearance_version` to be there.

## 5. Client repaints counts on an appearance update

`updateNodeAppearance` (`analysis_nodes.js:1064`) redraws the node from `node_data` but leaves the count
strip alone. Extract the count-painting block of `updateDirtyNode`'s `asyncUpdateNode` (the
`if (data.valid) { ... counts ... markLiveDataCount ... }` branch, `analysis_nodes.js:627-640`) into
`paintNodeCounts(node, data)` and use it from both places:

```js
function updateNodeAppearance(data) {
    const nodeId = data.attributes.node_id;
    const node = getNode(nodeId);
    updateNodeFromData(node, data);
    setupConnections(node);
    if (!node.attr("loading")) {
        messagePoller.observe_node(nodeId, "count", function(status) { paintNodeCounts(node, status); });
    }
}
```

`nodes_status` already returns every `NodeCount` label for the node's current version, so the repaint
needs no new endpoint. A node mid-reload keeps its spinner and gets its counts through the normal
dirty path.

Edit in `variantgrid/static_files/default_static/js/analysis_nodes.js`.

## 6. Tests

`analysis/tests/test_tagged_node_count.py`, built on `create_fake_variants` + an `AllVariantsNode`
(see `analysis/tests/test_variant_tags.py` for the analysis/tag scaffolding), with
`analysis.set_node_count_types([TOTAL, TAGGED])`:

- `get_extra_filters_q(analysis, TAGGED)` counts distinct variants: two tags on one variant plus one tag
  on another → 2. A tag on the same variant in a different analysis is excluded.
- Overlay, node loaded with a `TAGGED` row: first tag on a variant in the node → row +1 and
  `appearance_version` +1; second tag on the same variant → unchanged; a tag on a variant outside the
  node → unchanged; deleting one of two tags → unchanged; deleting the last → row −1.
- Tagged absent from the analysis's node count types → tagging changes neither `NodeCount` nor
  `appearance_version`.
- `_get_cached_label_count(TAGGED)` returns `None` for a node with `cloned_from` set, while another
  label still returns the clone's count.

Run `python3 manage.py test --keepdb analysis.tests.test_tagged_node_count analysis.tests.test_variant_tags
analysis.tests.test_node_count_provenance analysis.tests.test_built_in_filter_omim analysis.tests.test_urls`
- the last three cover the `get_extra_filters_q` signature change through the built-in filter node and
the grid.

## Files

| File | Change |
|---|---|
| `snpdb/models/models_enums.py` | `TAGGED`, `CHOICES`, `COLORS` |
| `snpdb/migrations/0224_*`, `analysis/migrations/0125_*` | generated choices migrations |
| `analysis/models/nodes/node_counts.py` | `get_extra_filters_q(analysis, extra_filters)` + Tagged branch |
| `analysis/models/nodes/filters/built_in_filter_node.py`, `analysis/grids.py` | new signature |
| `analysis/models/nodes/analysis_node.py` | skip `cloned_from` shortcut for Tagged |
| `analysis/variant_tag_node_counts.py` | `update_tagged_node_counts` |
| `analysis/signals/signal_handlers.py` | call it on tag create / delete |
| `variantgrid/static_files/default_static/js/analysis_nodes.js` | `paintNodeCounts`, repaint on appearance update |
| `analysis/tests/test_tagged_node_count.py` | tests above |
