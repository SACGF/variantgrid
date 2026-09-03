# Fold NodeCount into NodeVersion

Written by Claude Fable 5 (claude-fable-5), 2026-09-03

Node counts are a per-load snapshot of a NodeVersion, but they live in their own table: every node
load bulk-upserts one row per configured label and cascade-deletes the old version's rows, and
`variant_ids` sits on every row while only meaning something on the TOTAL row. Folding them into
NodeVersion leaves one row per version, cuts the churn on a hot create/delete path, and gives the
TagNode editor a place to snapshot its tag picker counts at load so the editor renders instantly
(issue #1820) - including in "All analyses" mode, whose count query is the slow one. Global-mode
counts going stale until the node is saved is exactly what the node's existing warning already
promises ("Global tags are a snapshot taken <date>... Press save to refresh").

## Data

```python
class NodeVersion(TimeStampedModel):
    """ Deleted once a node updates - all version-specific caches cascade from this """
    node = models.ForeignKey(AnalysisNode, on_delete=CASCADE)
    version = models.IntegerField(null=False)
    # {source_key: data_version} of the mutable tables this node read at load. Empty = deterministic
    live_data_sources = models.JSONField(default=dict)
    # The exact PK set stored at load, only for nodes <= ANALYSIS_NODE_STORE_ID_SIZE_MAX.
    # When present, load_data["counts"]["T"] == len(variant_ids)
    variant_ids = ArrayField(models.IntegerField(), null=True)
    # Products of the node's load:
    #   "counts":     {node count label: count} - the DAG badge counts, eg {"T": 1234, "C": 4, "tag_artefact": 3}.
    #                 Labels only ever live under this key, so they can never collide with the keys beside it
    #   "tag_counts": {tag: count} - TagNode only: the editor's tag picker, counted over the node's input
    load_data = models.JSONField(default=dict)

    class Meta:
        unique_together = ("node", "version")
```

`NodeCount` is deleted. `modified` (from TimeStampedModel) becomes the client's "a recount landed"
signal - after load, the only thing that ever touches the row is a tag recount, so it means exactly
what `NodeCount.modified` meant. Every count write goes through `.update()` / raw SQL, which bypass
`auto_now`, so each write sets `modified` explicitly - the current code already does this.

## Writers

- **`AnalysisNode.load`** (`analysis_node.py` ~1123, the `node_counts` step): replace the NodeCount
  `bulk_create(update_conflicts=True)` with one
  `NodeVersion.objects.filter(pk=...).update(load_data=..., variant_ids=..., modified=now())`.
  A full overwrite is correct here for the same reason the current comment gives: counts are a cache
  of a query against an immutable node_version, so a re-load can safely overwrite. The
  `live_data_sources` write earlier in `load` stays a separate targeted `.update()` - it must land
  before counting, and must not clobber a concurrent recount's counts.
  `load_data["tag_counts"]` is supplied by the node subclass: the base contributes nothing, TagNode
  contributes its picker counts (the existing `get_tag_counts()` body, which already counts the
  input scope and keeps configured tags at zero count).
- **`update_analysis_tag_node_counts`** (`node_utils.py`): the recount runs against versions the
  nodes are already on, concurrently with loads, so it must merge into `counts` without clobbering
  the labels and `tag_counts` it didn't compute - a DB-side deep merge, one UPDATE per node
  (an analysis has tens of nodes; the current version also builds each node's queryset in a loop):

  ```sql
  UPDATE analysis_nodeversion
  SET load_data = jsonb_set(load_data, '{counts}',
                            COALESCE(load_data->'counts', '{}'::jsonb) || %s::jsonb),
      modified = %s
  WHERE id = %s
  ```

  A node that bumped its version mid-count has had this row deleted, so the UPDATE touches 0 rows
  and the recount skips it naturally - the `IntegrityError` catch around today's bulk upsert goes.

## Readers

- **`nodes_status`** (`views_json.py`): collapses to a single NodeVersion query - `load_data`
  ("counts"), `modified` (as `counts_modified`) and `live_data_sources` all come off the one row,
  replacing the NodeCount join and the max-modified dedupe loop. TOTAL is still overridden from
  `node.count` as now.
- **`get_extra_filters_count`** (`node_counts.py`): `node_version.nodecount_set.filter(label=...)`
  becomes `node.node_version.load_data.get("counts", {}).get(extra_filters)`.
- **`get_cached_node_pks`** (`analysis_node.py`): reads `node.node_version.variant_ids`; the
  issue #546 explicit-PK substitution in `get_small_parent_arg_q_dict` is otherwise unchanged.
- **`BuiltInFilterNode._get_cached_label_count`**: the parent's count comes from the parent's
  `node_version.load_data` instead of `NodeCount.load_for_node`.
- **`analysis_input_stats`** (`views_analysis_settings.py`): the wall-clock span over
  `NodeCount.created` becomes the span over `NodeVersion.modified` for the analysis - a later tag
  recount skews it slightly, acceptable for a debug page.
- **TagNode editor** (`TagNodeView._get_tag_counts_context`): reads
  `node.node_version.load_data.get("tag_counts")` and sorts/labels at render (user tag sort order,
  `TagFilter.label`). When the key is absent - a version loaded before this deploys, or mid-reload -
  fall back to computing `get_tag_counts()` live, which is today's behaviour. `get_tag_counts()`
  itself moves to being called at load time.

`AnalysisNode.load_for_node` / `load_for_node_version` callers and the
`NodeCount.DoesNotExist` fallbacks they carry all reduce to dict lookups with defaults.

## Migration

One new migration on top (nothing here reshapes pushed migrations):

1. Add `variant_ids` and `load_data` to NodeVersion.
2. `RunPython`: fold existing NodeCount rows in, grouped by node_version - labels into
   `load_data["counts"]`, the TOTAL row's `variant_ids` into `NodeVersion.variant_ids`. Iterate in
   batches; deployments have a row per (live node version x configured label). Open analyses keep
   their badges without a forced reload.
3. Delete the NodeCount model.

Runs entirely under `manage.py migrate` - no ManualOperation needed.

## Tests

Update the NodeCount touchpoints to read/write NodeVersion instead, keeping the behaviours they
cover: `test_tag_node_counts` (recount values, `counts_modified` signal - now `NodeVersion.modified`,
editor pills), `test_explicit_pk_substitution`, `test_grid_known_count`, `test_node_grid_auto_load`,
`test_node_count_provenance`, `test_clone_nodes`. Add one test that a tag recount merging into
`counts` preserves an existing `tag_counts` key.

## Finishing

- `scripts/vg map` (model changed) and update the NodeCount mentions in `analysis/CLAUDE.md`
  ("Small parents are inlined ... from the pks stored in NodeCount", tag recount gotcha) and
  `__analysis_readme.md`.
