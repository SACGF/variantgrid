# #1838 — TagNode load: take tag counting out of the load, and time the phases

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-08
Status: approved

[#1838](https://github.com/SACGF/variantgrid/issues/1838): a prod TagNode recorded `load_seconds = 149.96` on an
analysis with **0 VariantTags**. The node's own filter is `Q(pk__in=[])`, which Django short-circuits without SQL,
so the time is in the rest of `AnalysisNode.load()` - and `load_seconds` is one number, so nothing says where.

Two changes:

1. **The tag picker stops being part of the load.** Its counts are hints beside the pills; they come from the
   tagging table alone and are worked out when the editor renders. `TagNode` no longer counts anything at load,
   snapshots nothing, and never builds a queryset over `snpdb_variant` to decorate a picker.
2. **Loads record per-phase timings**, so the next slow load on prod arrives with its phase named.

## 1. What was confirmed on vg-test2 (200k `analysis_varianttag` rows, 28 TagNodes)

| node | shape | recorded `load_seconds` | re-timed now |
|---|---|---|---|
| 225, 254, 221, 252 | hidden "Tagged Variants" source, THIS_ANALYSIS, 0 tags in analysis | 5.3 - 7.3 s | `get_tag_counts` 0.016 s, label counts 0.014 s |
| 404, 490 | global source (admin) | 0.3 - 2.6 s | `get_tag_counts` 0.25 s, 25 queries |
| 488 | global source (non-admin, guardian subquery) | 0.6 s | `get_tag_counts` 1.06 s, 26 queries |

Nothing in the load path takes seconds when re-run warm, yet the empty hidden nodes recorded 5-7 s. They are
created by `analysis/models/nodes/filters/tag_node.py:TagNode.get_analysis_tags_node`, which runs `update_node_task`
synchronously inside the view (`analysis/views/views.py` line 167) - cold cache, a busy box or a lock wait inside
`node_counts()` would all land in `load_seconds` and look the same. The prod 150 s is the same blind spot at SA Path
scale (412k taggings). Hence §4.

What the current picker costs, per TagNode load, whatever the analysis holds:

- The tag list `Tag.objects.filter(Q(varianttag__analysis=...) | Q(tagnodetag__tag_node=...)).distinct()` plans as
  a seq scan of the whole tagging table (29 ms cold here; the `EXPLAIN ANALYZE` is in the issue thread).
- In `TagNodeMode.ALL_TAGS`, one `COUNT(v.id) FILTER (WHERE va.allele_id IN (<all taggings>))` per live tag over
  `snpdb_variant ⟕ snpdb_variantallele`, and for a source node the outer scope is every variant in the build.

For comparison, `VariantTag.objects.filter(unresolved).values("tag_id").annotate(Count("id"))` over all 200k rows
is 33 ms (parallel seq scan, one hash aggregate), and the same restricted to one analysis is an index scan.

## 2. Data

No model changes. `NodeVersion.load_data` (JSON, written by `node_counts`) loses `tag_counts` and gains `timings`:

```python
load_data = {
    "counts": {label: count},           # existing
    "timings": {                        # new - seconds per load phase, rounded to 3 dp
        "load": float,                  #   _load()
        "live_data_sources": float,     #   get_live_data_sources()
        "counts": float,                #   label counts (get_node_counts_and_labels_dict)
        "variant_ids": float,           #   _get_variant_ids_to_store()
        "load_data": float,             #   _get_load_data()
    },
}
```

Existing NodeVersion rows keep their `tag_counts` key; nothing reads it after this change, and the rows are dropped
by `delete_analysis_old_node_versions` as nodes reload. No migration.

## 3. The picker: counts from the tagging table, at render time

`TagNode.get_tag_counts()` becomes one grouped query over `VariantTag`, no queryset, no parent, no `snpdb_variant`:

```python
def get_tag_counts(self) -> dict[str, int]:
    """ {tag: taggings} for the editor's tag picker - a hint beside each pill, not the node's count.
        Local mode counts this analysis's taggings; global mode counts every tagging the user can see,
        whatever the analysis or build - the node's own filter decides what gets through """
    if self.mode == TagNodeMode.ALL_TAGS:
        tags_qs = VariantTag.filter_for_user(self.analysis.user)
    else:
        tags_qs = VariantTag.objects.filter(analysis=self.analysis)
    if not self.include_resolved:
        tags_qs = tags_qs.filter(VariantTag.unresolved_q())
    tag_counts = dict(tags_qs.values_list("tag_id").annotate(n=Count("id")).values_list("tag_id", "n"))
    # A configured tag always keeps its pill - dropping it would silently drop the tag on the next save
    for tag_id in self.tag_ids:
        tag_counts.setdefault(tag_id, 0)
    if self.mode == TagNodeMode.ALL_TAGS:
        for tag_id in Tag.objects.filter(retired__isnull=True).values_list("pk", flat=True):
            tag_counts.setdefault(tag_id, 0)  # every live tag is pickable in global mode
    return tag_counts
```

What changes in meaning: the number beside a pill is **taggings in scope**, not "variants in this node's input
carrying the tag". In local mode those are the same thing in practice (one tagging per variant per sample in an
analysis, `varianttag_one_per_sample_in_analysis`). In global mode it is system-wide and build-agnostic - the same
figure the analyses list and the variant tags page already show through `tag_counts_summary`. The
`tagged_within_days` cutoff stays out of it, as now.

Consequences, all deletions:

- `TagNode._get_load_data` goes, and with it the `tag_counts` key in `load_data` and the comment on it in
  `analysis/models/nodes/analysis_node.py` (line ~1443). The base `_get_load_data` hook stays for other nodes.
- `analysis/views/nodes/node_views.py:TagNodeView._get_tag_counts_context` calls `node.get_tag_counts()` directly.
  The `NodeStatus.is_ready` gate and the `show_tag_counts` flag go: the picker no longer depends on the load, so the
  editor template's "The tag picker appears once the node has loaded" branch goes too.
- The `NonFatalNodeError` catch in `get_tag_counts` goes - there is no parent queryset to fail on.
- The `TagNodeMode.ALL_TAGS` snapshot warning in `TagNode.get_warnings` stays: it is about the node's variants, which
  are still a snapshot of the taggings at load. The picker is now live, which is fine for a hint.
- `analysis/CLAUDE.md`: the line "TagNode snapshots its editor's `tag_counts` there, as counting them in global mode
  is slow" comes out of the `_load`/`load_data` bullet.

No global count cache is needed: the query is one aggregate over a table indexed on `analysis_id` and small enough
to hash-aggregate in tens of milliseconds at 200k rows, and it runs on editor open, not on every load.

## 4. Phase timing in `AnalysisNode.load()` / `node_counts()`

- `load()` times `_load()`; `node_counts()` times each of its steps with `time.perf_counter()` and puts the dict
  into `load_data["timings"]` before the single `NodeVersion.update(load_data=...)` it already does. No extra write.
- One `logging.warning("Node %d.%d slow load %.1fs: %s", ...)` with the timings dict when the total exceeds
  `settings.ANALYSIS_NODE_SLOW_LOAD_SECONDS` (new setting, default 30). Rollbar picks warnings up, so the next 150 s
  arrives with its phase breakdown attached.
- `manage.py profile_analysis_nodes` already prints `cached_load_seconds`; add the `timings` dict to its per-node
  row so `--analysis <id>` on prod answers the issue's "to confirm on prod" step without `--rerun`.

## 5. Tests (`analysis/tests/test_tag_node_counts.py`, `analysis/tests/test_variant_tags.py`)

Keep, unchanged - they assert the picker's contents and still hold under the new definition:

- `TestTagNodeEditorCounts` (pill per tag in the analysis, configured tag keeps its pill at zero, other tags stay
  pickable, count of 1 rendered) - `test_counts_cover_the_input_so_other_tags_stay_pickable` keeps its assertion;
  its docstring changes to say the picker counts the scope's taggings.
- `TestGlobalTagNodeCounts.test_counts_a_tag_made_in_another_analysis`, `test_ignores_the_tagged_within_days_cutoff`,
  `test_local_mode_does_not_count_another_analysis_tag` - the global-mode ones now also see `{other_tag: 0}` pinned
  for every live tag, so the assertions become `assertEqual(1, node.get_tag_counts()[self.tag.pk])`.
- `test_variant_tags.py` line 337: `include_resolved` counts the resolved tagging too - unchanged.

Change: `TestTagNodeCountValues.test_recount_keeps_the_rest_of_the_load_data` (line ~148) seeds `tag_counts` as the
"rest" - swap the key for `timings` so it still proves the jsonb merge keeps the other keys.

Add:

- `TestTagNodeEditorCounts`: global mode over an empty analysis renders a pill for a live tag tagged in another
  analysis and issues no query naming `snpdb_variant` (`CaptureQueriesContext`).
- One `load()` test asserting `node_version.load_data["timings"]` has the five keys - the phase boundaries are our
  logic.

Drop nothing else.

## 6. Files

- `analysis/models/nodes/filters/tag_node.py` - `get_tag_counts` rewritten, `_get_load_data` removed (§3)
- `analysis/views/nodes/node_views.py` - `_get_tag_counts_context` (§3)
- `analysis/templates/analysis/node_editors/tagnode_editor.html` - drop the `show_tag_counts` branches (§3)
- `analysis/models/nodes/analysis_node.py` - `load`, `node_counts`, `load_data` docstring (§4)
- `analysis/management/commands/profile_analysis_nodes.py` - emit `timings` (§4)
- `variantgrid/settings/components/default_settings.py` - `ANALYSIS_NODE_SLOW_LOAD_SECONDS`
- `analysis/CLAUDE.md` - amend the `load_data` bullet; add under Gotchas that the picker counts taggings from
  `analysis_varianttag` at render, and that `load_data["timings"]` is where to look when `load_seconds` is big
- tests as §5
- `scripts/vg map` after the setting is added

## 7. Verification

- `python3 manage.py test --keepdb analysis.tests.test_tag_node_counts analysis.tests.test_variant_tags`
- `vg page node_view --kwargs ... --queries` on the editors of 225 (empty local source), 404/490 (global source,
  admin) and 488 (global source, guardian-scoped user): no query names `snpdb_variant`; pills and counts render
- Reload 404 and 490 via their editors on test.variantgrid.com: `load_data` has `timings`, no `tag_counts`;
  `load_seconds` drops from 0.3-2.6 s to the label-count cost alone
- `manage.py profile_analysis_nodes --analysis 39` shows the timings column

## 8. Open question for prod

The prod node id and mode are still unknown. Once §4 lands and the node is reloaded, `load_data["timings"]` on its
NodeVersion says which phase. With the picker out of the load, a TagNode's load is `node_counts()` like any other
node's, so a repeat points at the label-count query (`get_node_counts_and_labels_dict`) or at contention on the
box, not at tags.
